/**
 * This file is part of ORB-SLAM3
 *
 * Copyright (C) 2017-2021 Carlos Campos, Richard Elvira, Juan J. Gómez Rodríguez, José M.M. Montiel and Juan D. Tardós, University of Zaragoza.
 * Copyright (C) 2014-2016 Raúl Mur-Artal, José M.M. Montiel and Juan D. Tardós, University of Zaragoza.
 *
 * ORB-SLAM3 is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ORB-SLAM3 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with ORB-SLAM3.
 * If not, see <http://www.gnu.org/licenses/>.
 */

#include "LocalMapping.h"
#include "Converter.h"
#include "GeometricTools.h"
#include "LoopClosing.h"
#include "ORBmatcher.h"
#include "Optimizer.h"

#include <chrono>
#include <mutex>

namespace ORB_SLAM3 {

LocalMapping::LocalMapping(System* pSys, Atlas* pAtlas, const float bMonocular, bool bInertial, const string& _strSeqName)
    : mpSystem(pSys)
    , mbMonocular(bMonocular)
    , mbInertial(bInertial)
    , mbResetRequested(false)
    , mbResetRequestedActiveMap(false)
    , mbFinishRequested(false)
    , mbFinished(true)
    , mpAtlas(pAtlas)
    , bInitializing(false)
    , mbAbortBA(false)
    , mbStopped(false)
    , mbStopRequested(false)
    , mbNotStop(false)
    , mbAcceptKeyFrames(true)
    , mIdxInit(0)
    , mScale(1.0)
    , mInitSect(0)
    , mbNotBA1(true)
    , mbNotBA2(true)
    , mIdxIteration(0)
    , infoInertial(Eigen::MatrixXd::Zero(9, 9))
{
    mnMatchesInliers = 0;

    mbBadImu = false;

    mTinit = 0.f;

    mNumLM = 0;
    mNumKFCulling = 0;

#ifdef REGISTER_TIMES
    nLBA_exec = 0;
    nLBA_abort = 0;
#endif
}

void LocalMapping::SetLoopCloser(LoopClosing* pLoopCloser)
{
    mpLoopCloser = pLoopCloser;
}

void LocalMapping::SetTracker(Tracking* pTracker)
{
    mpTracker = pTracker;
}

/**
 * @brief 局部地图线程主函数
 */
void LocalMapping::Run()
{
    // 标记状态，表示当前run函数正在运行，尚未结束
    mbFinished = false;

    while (1) {
        // Tracking will see that Local Mapping is busy
        // Step 1 告诉Tracking，LocalMapping正处于繁忙状态，请不要给我发送关键帧打扰我
        SetAcceptKeyFrames(false);

        // Check if there are keyframes in the queue
        // 等待处理的关键帧列表不为空 并且imu正常
        if (CheckNewKeyFrames() && !mbBadImu) {
#ifdef REGISTER_TIMES
            double timeLBA_ms = 0;
            double timeKFCulling_ms = 0;

            std::chrono::steady_clock::time_point time_StartProcessKF = std::chrono::steady_clock::now();
#endif
            // BoW conversion and insertion in Map
            // Step 2 处理列表中的关键帧，包括计算BoW、更新观测、描述子、共视图，插入到地图等
            ProcessNewKeyFrame();
#ifdef REGISTER_TIMES
            std::chrono::steady_clock::time_point time_EndProcessKF = std::chrono::steady_clock::now();

            double timeProcessKF = std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(time_EndProcessKF - time_StartProcessKF).count();
            vdKFInsert_ms.push_back(timeProcessKF);
#endif

            // Check recent MapPoints
            // Step 3 根据地图点的观测情况剔除质量不好的地图点
            MapPointCulling();
#ifdef REGISTER_TIMES
            std::chrono::steady_clock::time_point time_EndMPCulling = std::chrono::steady_clock::now();

            double timeMPCulling = std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(time_EndMPCulling - time_EndProcessKF).count();
            vdMPCulling_ms.push_back(timeMPCulling);
#endif

            // Triangulate new MapPoints
            // Step 4 当前关键帧与相邻关键帧通过三角化产生新的地图点，使得跟踪更稳
            CreateNewMapPoints();

            mbAbortBA = false;

            // 已经处理完队列中的最后的一个关键帧
            if (!CheckNewKeyFrames()) {
                // Find more matches in neighbor keyframes and fuse point duplications
                //  Step 5 检查并融合当前关键帧与相邻关键帧帧（两级相邻）中重复的地图点
                // 先完成相邻关键帧与当前关键帧的地图点的融合（在相邻关键帧中查找当前关键帧的地图点），
                // 再完成当前关键帧与相邻关键
                SearchInNeighbors();
            }

#ifdef REGISTER_TIMES
            std::chrono::steady_clock::time_point time_EndMPCreation = std::chrono::steady_clock::now();

            double timeMPCreation = std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(time_EndMPCreation - time_EndMPCulling).count();
            vdMPCreation_ms.push_back(timeMPCreation);
#endif

            bool b_doneLBA = false;
            int num_FixedKF_BA = 0;
            int num_OptKF_BA = 0;
            int num_MPs_BA = 0;
            int num_edges_BA = 0;

            // 已经处理完队列中的最后的一个关键帧，并且闭环检测没有请求停止LocalMapping
            if (!CheckNewKeyFrames() && !stopRequested()) {
                // 当前地图中关键帧数目大于2个
                if (mpAtlas->KeyFramesInMap() > 2) {
                    // Step 6.1 处于IMU模式并且当前关键帧所在的地图已经完成IMU初始化
                    if (mbInertial && mpCurrentKeyFrame->GetMap()->isImuInitialized()) {
                        // 计算上一关键帧到当前关键帧相机光心的距离 + 上上关键帧到上一关键帧相机光心的距离
                        float dist = (mpCurrentKeyFrame->mPrevKF->GetCameraCenter() - mpCurrentKeyFrame->GetCameraCenter()).norm() + (mpCurrentKeyFrame->mPrevKF->mPrevKF->GetCameraCenter() - mpCurrentKeyFrame->mPrevKF->GetCameraCenter()).norm();

                        // 如果距离大于5厘米，记录当前KF和上一KF时间戳的差，累加到mTinit
                        if (dist > 0.05)
                            mTinit += mpCurrentKeyFrame->mTimeStamp - mpCurrentKeyFrame->mPrevKF->mTimeStamp;
                        // 当前关键帧所在的地图尚未完成IMU BA2（IMU第三阶段初始化）
                        if (!mpCurrentKeyFrame->GetMap()->GetIniertialBA2()) {
                            // 如果累计时间差小于10s 并且 距离小于2厘米，认为运动幅度太小，不足以初始化IMU，将mbBadImu设置为true
                            if ((mTinit < 10.f) && (dist < 0.02)) {
                                cout << "Not enough motion for initializing. Reseting..." << endl;
                                unique_lock<mutex> lock(mMutexReset);
                                mbResetRequestedActiveMap = true;
                                mpMapToReset = mpCurrentKeyFrame->GetMap();
                                mbBadImu = true; // 在跟踪线程里会重置当前活跃地图
                            }
                        }

                        // 判断成功跟踪匹配的点数是否足够多
                        // 条件---------1.1、跟踪成功的内点数目大于75-----1.2、并且是单目--或--2.1、跟踪成功的内点数目大于100-----2.2、并且不是单目
                        bool bLarge = ((mpTracker->GetMatchesInliers() > 75) && mbMonocular) || ((mpTracker->GetMatchesInliers() > 100) && !mbMonocular);
                        // 局部地图+IMU一起优化，优化关键帧位姿、地图点、IMU参数
                        // 注意这里的第二个参数是按地址传递的,当这里的 mbAbortBA 状态发生变化时，能够及时执行/停止BA
                        Optimizer::LocalInertialBA(mpCurrentKeyFrame, &mbAbortBA, mpCurrentKeyFrame->GetMap(), num_FixedKF_BA, num_OptKF_BA, num_MPs_BA, num_edges_BA, bLarge, !mpCurrentKeyFrame->GetMap()->GetIniertialBA2());
                        b_doneLBA = true;
                    } else {
                        // 注意这里的第二个参数是按地址传递的,当这里的 mbAbortBA 状态发生变化时，能够及时执行/停止BA
                        Optimizer::LocalBundleAdjustment(mpCurrentKeyFrame, &mbAbortBA, mpCurrentKeyFrame->GetMap(), num_FixedKF_BA, num_OptKF_BA, num_MPs_BA, num_edges_BA);
                        b_doneLBA = true;
                    }
                }
#ifdef REGISTER_TIMES
                std::chrono::steady_clock::time_point time_EndLBA = std::chrono::steady_clock::now();

                if (b_doneLBA) {
                    timeLBA_ms = std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(time_EndLBA - time_EndMPCreation).count();
                    vdLBA_ms.push_back(timeLBA_ms);

                    nLBA_exec += 1;
                    if (mbAbortBA) {
                        nLBA_abort += 1;
                    }
                    vnLBA_edges.push_back(num_edges_BA);
                    vnLBA_KFopt.push_back(num_OptKF_BA);
                    vnLBA_KFfixed.push_back(num_FixedKF_BA);
                    vnLBA_MPs.push_back(num_MPs_BA);
                }

#endif
                // Initialize IMU here
                // Step 7 初始化IMU，确定scale，gravity，bias
                if (!mpCurrentKeyFrame->GetMap()->isImuInitialized() && mbInertial) {
                    if (mbMonocular)
                        InitializeIMU(1e2, 1e10, true);
                    else
                        InitializeIMU(1e2, 1e5, true);
                }

                // Check redundant local Keyframes
                // Step 8 检测并剔除当前帧相邻的关键帧中冗余的关键帧
                // 冗余的判定：该关键帧的90%的地图点可以被其它关键帧观测到
                KeyFrameCulling();

#ifdef REGISTER_TIMES
                std::chrono::steady_clock::time_point time_EndKFCulling = std::chrono::steady_clock::now();

                timeKFCulling_ms = std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(time_EndKFCulling - time_EndLBA).count();
                vdKFCulling_ms.push_back(timeKFCulling_ms);
#endif
                // Step 9: IMU参数分时间进行优化求解
                if ((mTinit < 50.0f) && mbInertial) {
                    if (mpCurrentKeyFrame->GetMap()->isImuInitialized() && mpTracker->mState == Tracking::OK) // Enter here everytime local-mapping is called
                    {
                        if (!mpCurrentKeyFrame->GetMap()->GetIniertialBA1()) {
                            if (mTinit > 5.0f) {
                                cout << "start VIBA 1" << endl;
                                mpCurrentKeyFrame->GetMap()->SetIniertialBA1();
                                if (mbMonocular)
                                    InitializeIMU(1.f, 1e5, true);
                                else
                                    InitializeIMU(1.f, 1e5, true);

                                cout << "end VIBA 1" << endl;
                            }
                        } else if (!mpCurrentKeyFrame->GetMap()->GetIniertialBA2()) {
                            if (mTinit > 15.0f) {
                                cout << "start VIBA 2" << endl;
                                mpCurrentKeyFrame->GetMap()->SetIniertialBA2();
                                if (mbMonocular)
                                    InitializeIMU(0.f, 0.f, true);
                                else
                                    InitializeIMU(0.f, 0.f, true);

                                cout << "end VIBA 2" << endl;
                            }
                        }

                        // scale refinement
                        if (((mpAtlas->KeyFramesInMap()) <= 200) && ((mTinit > 25.0f && mTinit < 25.5f) || (mTinit > 35.0f && mTinit < 35.5f) || (mTinit > 45.0f && mTinit < 45.5f) || (mTinit > 55.0f && mTinit < 55.5f) || (mTinit > 65.0f && mTinit < 65.5f) || (mTinit > 75.0f && mTinit < 75.5f))) {
                            if (mbMonocular)
                                ScaleRefinement();
                        }
                    }
                }
            }

#ifdef REGISTER_TIMES
            vdLBASync_ms.push_back(timeKFCulling_ms);
            vdKFCullingSync_ms.push_back(timeKFCulling_ms);
#endif

            mpLoopCloser->InsertKeyFrame(mpCurrentKeyFrame);

#ifdef REGISTER_TIMES
            std::chrono::steady_clock::time_point time_EndLocalMap = std::chrono::steady_clock::now();

            double timeLocalMap = std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(time_EndLocalMap - time_StartProcessKF).count();
            vdLMTotal_ms.push_back(timeLocalMap);
#endif
        } else if (Stop() && !mbBadImu) {
            // Safe area to stop
            while (isStopped() && !CheckFinish()) {
                usleep(3000);
            }
            if (CheckFinish())
                break;
        }

        ResetIfRequested();

        // Tracking will see that Local Mapping is busy
        SetAcceptKeyFrames(true);

        if (CheckFinish())
            break;

        usleep(3000);
    }

    SetFinish();
}
/**
 * @brief 插入关键帧,由外部（Tracking）线程调用;这里只是插入到列表中,等待线程主函数对其进行处理
 * @param pKF 新的关键帧
 */
void LocalMapping::InsertKeyFrame(KeyFrame* pKF)
{
    unique_lock<mutex> lock(mMutexNewKFs);
    mlNewKeyFrames.push_back(pKF);
    mbAbortBA = true;
}

bool LocalMapping::CheckNewKeyFrames()
{
    unique_lock<mutex> lock(mMutexNewKFs);
    return (!mlNewKeyFrames.empty());
}
/**
 * @brief 处理列表中的关键帧，包括计算BoW、更新观测、描述子、共视图，插入到地图等
 */
void LocalMapping::ProcessNewKeyFrame()
{
    // Step 1：从缓冲队列中取出一帧关键帧
    // 该关键帧队列是Tracking线程向LocalMapping中插入的关键帧组成
    {
        unique_lock<mutex> lock(mMutexNewKFs);
        mpCurrentKeyFrame = mlNewKeyFrames.front();
        mlNewKeyFrames.pop_front();
    }

    // Compute Bags of Words structures
    mpCurrentKeyFrame->ComputeBoW();

    // Associate MapPoints to the new keyframe and update normal and descriptor
    // Step 3：当前处理关键帧中有效的地图点，更新normal，描述子等信息
    // TrackLocalMap中和当前帧新匹配上的地图点和当前关键帧进行关联绑定
    const vector<MapPoint*> vpMapPointMatches = mpCurrentKeyFrame->GetMapPointMatches();

    // 对当前处理的这个关键帧中的所有的地图点展开遍历
    for (size_t i = 0; i < vpMapPointMatches.size(); i++) {
        MapPoint* pMP = vpMapPointMatches[i];
        if (pMP) {
            if (!pMP->isBad()) {
                if (!pMP->IsInKeyFrame(mpCurrentKeyFrame)) {
                    // 如果地图点不是来自当前帧的观测，为当前地图点添加观测
                    pMP->AddObservation(mpCurrentKeyFrame, i);
                    // 获得该点的平均观测方向和观测距离范围
                    pMP->UpdateNormalAndDepth();
                    // 更新地图点的最佳描述子
                    pMP->ComputeDistinctiveDescriptors();
                } else // this can only happen for new stereo points inserted by the Tracking
                {
                    // 如果当前帧中已经包含了这个地图点,但是这个地图点中却没有包含这个关键帧的信息
                    // 这些地图点可能来自双目或RGBD跟踪过程中新生成的地图点，或者是CreateNewMapPoints 中通过三角化产生
                    // 将上述地图点放入mlpRecentAddedMapPoints，等待后续MapPointCulling函数的检验
                    mlpRecentAddedMapPoints.push_back(pMP);
                }
            }
        }
    }

    // Update links in the Covisibility Graph
    mpCurrentKeyFrame->UpdateConnections();

    // Insert Keyframe in Map
    mpAtlas->AddKeyFrame(mpCurrentKeyFrame);
}

void LocalMapping::EmptyQueue()
{
    while (CheckNewKeyFrames())
        ProcessNewKeyFrame();
}

void LocalMapping::MapPointCulling()
{
    // Check Recent Added MapPoints
    list<MapPoint*>::iterator lit = mlpRecentAddedMapPoints.begin();
    const unsigned long int nCurrentKFid = mpCurrentKeyFrame->mnId;

    int nThObs;
    if (mbMonocular)
        nThObs = 2;
    else
        nThObs = 3;
    const int cnThObs = nThObs;

    int borrar = mlpRecentAddedMapPoints.size();

    while (lit != mlpRecentAddedMapPoints.end()) {
        MapPoint* pMP = *lit;

        if (pMP->isBad())
            lit = mlpRecentAddedMapPoints.erase(lit);
        else if (pMP->GetFoundRatio() < 0.25f) {
            pMP->SetBadFlag();
            lit = mlpRecentAddedMapPoints.erase(lit);
        } else if (((int)nCurrentKFid - (int)pMP->mnFirstKFid) >= 2 && pMP->Observations() <= cnThObs) {
            pMP->SetBadFlag();
            lit = mlpRecentAddedMapPoints.erase(lit);
        } else if (((int)nCurrentKFid - (int)pMP->mnFirstKFid) >= 3)
            lit = mlpRecentAddedMapPoints.erase(lit);
        else {
            lit++;
            borrar--;
        }
    }
}

void LocalMapping::CreateNewMapPoints()
{
    // Retrieve neighbor keyframes in covisibility graph
    int nn = 10;
    // For stereo inertial case
    if (mbMonocular)
        nn = 30;
    vector<KeyFrame*> vpNeighKFs = mpCurrentKeyFrame->GetBestCovisibilityKeyFrames(nn);

    if (mbInertial) {
        KeyFrame* pKF = mpCurrentKeyFrame;
        int count = 0;
        while ((vpNeighKFs.size() <= nn) && (pKF->mPrevKF) && (count++ < nn)) {
            vector<KeyFrame*>::iterator it = std::find(vpNeighKFs.begin(), vpNeighKFs.end(), pKF->mPrevKF);
            if (it == vpNeighKFs.end())
                vpNeighKFs.push_back(pKF->mPrevKF);
            pKF = pKF->mPrevKF;
        }
    }

    float th = 0.6f;

    ORBmatcher matcher(th, false);

    Sophus::SE3<float> sophTcw1 = mpCurrentKeyFrame->GetPose();
    Eigen::Matrix<float, 3, 4> eigTcw1 = sophTcw1.matrix3x4();
    Eigen::Matrix<float, 3, 3> Rcw1 = eigTcw1.block<3, 3>(0, 0);
    Eigen::Matrix<float, 3, 3> Rwc1 = Rcw1.transpose();
    Eigen::Vector3f tcw1 = sophTcw1.translation();
    Eigen::Vector3f Ow1 = mpCurrentKeyFrame->GetCameraCenter();

    const float& fx1 = mpCurrentKeyFrame->fx;
    const float& fy1 = mpCurrentKeyFrame->fy;
    const float& cx1 = mpCurrentKeyFrame->cx;
    const float& cy1 = mpCurrentKeyFrame->cy;
    const float& invfx1 = mpCurrentKeyFrame->invfx;
    const float& invfy1 = mpCurrentKeyFrame->invfy;

    const float ratioFactor = 1.5f * mpCurrentKeyFrame->mfScaleFactor;
    int countStereo = 0;
    int countStereoGoodProj = 0;
    int countStereoAttempt = 0;
    int totalStereoPts = 0;
    // Search matches with epipolar restriction and triangulate
    for (size_t i = 0; i < vpNeighKFs.size(); i++) {
        if (i > 0 && CheckNewKeyFrames())
            return;

        KeyFrame* pKF2 = vpNeighKFs[i];

        GeometricCamera *pCamera1 = mpCurrentKeyFrame->mpCamera, *pCamera2 = pKF2->mpCamera;

        // Check first that baseline is not too short
        Eigen::Vector3f Ow2 = pKF2->GetCameraCenter();
        Eigen::Vector3f vBaseline = Ow2 - Ow1;
        const float baseline = vBaseline.norm();

        if (!mbMonocular) {
            if (baseline < pKF2->mb)
                continue;
        } else {
            const float medianDepthKF2 = pKF2->ComputeSceneMedianDepth(2);
            const float ratioBaselineDepth = baseline / medianDepthKF2;

            if (ratioBaselineDepth < 0.01)
                continue;
        }

        // Search matches that fullfil epipolar constraint
        vector<pair<size_t, size_t>> vMatchedIndices;
        bool bCoarse = mbInertial && mpTracker->mState == Tracking::RECENTLY_LOST && mpCurrentKeyFrame->GetMap()->GetIniertialBA2();

        matcher.SearchForTriangulation(mpCurrentKeyFrame, pKF2, vMatchedIndices, false, bCoarse);

        Sophus::SE3<float> sophTcw2 = pKF2->GetPose();
        Eigen::Matrix<float, 3, 4> eigTcw2 = sophTcw2.matrix3x4();
        Eigen::Matrix<float, 3, 3> Rcw2 = eigTcw2.block<3, 3>(0, 0);
        Eigen::Matrix<float, 3, 3> Rwc2 = Rcw2.transpose();
        Eigen::Vector3f tcw2 = sophTcw2.translation();

        const float& fx2 = pKF2->fx;
        const float& fy2 = pKF2->fy;
        const float& cx2 = pKF2->cx;
        const float& cy2 = pKF2->cy;
        const float& invfx2 = pKF2->invfx;
        const float& invfy2 = pKF2->invfy;

        // Triangulate each match
        const int nmatches = vMatchedIndices.size();
        for (int ikp = 0; ikp < nmatches; ikp++) {
            const int& idx1 = vMatchedIndices[ikp].first;
            const int& idx2 = vMatchedIndices[ikp].second;

            const cv::KeyPoint& kp1 = (mpCurrentKeyFrame->NLeft == -1) ? mpCurrentKeyFrame->mvKeysUn[idx1]
                : (idx1 < mpCurrentKeyFrame->NLeft)                    ? mpCurrentKeyFrame->mvKeys[idx1]
                                                                       : mpCurrentKeyFrame->mvKeysRight[idx1 - mpCurrentKeyFrame->NLeft];
            const float kp1_ur = mpCurrentKeyFrame->mvuRight[idx1];
            bool bStereo1 = (!mpCurrentKeyFrame->mpCamera2 && kp1_ur >= 0);
            const bool bRight1 = (mpCurrentKeyFrame->NLeft == -1 || idx1 < mpCurrentKeyFrame->NLeft) ? false
                                                                                                     : true;

            const cv::KeyPoint& kp2 = (pKF2->NLeft == -1) ? pKF2->mvKeysUn[idx2]
                : (idx2 < pKF2->NLeft)                    ? pKF2->mvKeys[idx2]
                                                          : pKF2->mvKeysRight[idx2 - pKF2->NLeft];

            const float kp2_ur = pKF2->mvuRight[idx2];
            bool bStereo2 = (!pKF2->mpCamera2 && kp2_ur >= 0);
            const bool bRight2 = (pKF2->NLeft == -1 || idx2 < pKF2->NLeft) ? false
                                                                           : true;

            if (mpCurrentKeyFrame->mpCamera2 && pKF2->mpCamera2) {
                if (bRight1 && bRight2) {
                    sophTcw1 = mpCurrentKeyFrame->GetRightPose();
                    Ow1 = mpCurrentKeyFrame->GetRightCameraCenter();

                    sophTcw2 = pKF2->GetRightPose();
                    Ow2 = pKF2->GetRightCameraCenter();

                    pCamera1 = mpCurrentKeyFrame->mpCamera2;
                    pCamera2 = pKF2->mpCamera2;
                } else if (bRight1 && !bRight2) {
                    sophTcw1 = mpCurrentKeyFrame->GetRightPose();
                    Ow1 = mpCurrentKeyFrame->GetRightCameraCenter();

                    sophTcw2 = pKF2->GetPose();
                    Ow2 = pKF2->GetCameraCenter();

                    pCamera1 = mpCurrentKeyFrame->mpCamera2;
                    pCamera2 = pKF2->mpCamera;
                } else if (!bRight1 && bRight2) {
                    sophTcw1 = mpCurrentKeyFrame->GetPose();
                    Ow1 = mpCurrentKeyFrame->GetCameraCenter();

                    sophTcw2 = pKF2->GetRightPose();
                    Ow2 = pKF2->GetRightCameraCenter();

                    pCamera1 = mpCurrentKeyFrame->mpCamera;
                    pCamera2 = pKF2->mpCamera2;
                } else {
                    sophTcw1 = mpCurrentKeyFrame->GetPose();
                    Ow1 = mpCurrentKeyFrame->GetCameraCenter();

                    sophTcw2 = pKF2->GetPose();
                    Ow2 = pKF2->GetCameraCenter();

                    pCamera1 = mpCurrentKeyFrame->mpCamera;
                    pCamera2 = pKF2->mpCamera;
                }
                eigTcw1 = sophTcw1.matrix3x4();
                Rcw1 = eigTcw1.block<3, 3>(0, 0);
                Rwc1 = Rcw1.transpose();
                tcw1 = sophTcw1.translation();

                eigTcw2 = sophTcw2.matrix3x4();
                Rcw2 = eigTcw2.block<3, 3>(0, 0);
                Rwc2 = Rcw2.transpose();
                tcw2 = sophTcw2.translation();
            }

            // Check parallax between rays
            Eigen::Vector3f xn1 = pCamera1->unprojectEig(kp1.pt);
            Eigen::Vector3f xn2 = pCamera2->unprojectEig(kp2.pt);

            Eigen::Vector3f ray1 = Rwc1 * xn1;
            Eigen::Vector3f ray2 = Rwc2 * xn2;
            const float cosParallaxRays = ray1.dot(ray2) / (ray1.norm() * ray2.norm());

            float cosParallaxStereo = cosParallaxRays + 1;
            float cosParallaxStereo1 = cosParallaxStereo;
            float cosParallaxStereo2 = cosParallaxStereo;

            if (bStereo1)
                cosParallaxStereo1 = cos(2 * atan2(mpCurrentKeyFrame->mb / 2, mpCurrentKeyFrame->mvDepth[idx1]));
            else if (bStereo2)
                cosParallaxStereo2 = cos(2 * atan2(pKF2->mb / 2, pKF2->mvDepth[idx2]));

            if (bStereo1 || bStereo2)
                totalStereoPts++;

            cosParallaxStereo = min(cosParallaxStereo1, cosParallaxStereo2);

            Eigen::Vector3f x3D;

            bool goodProj = false;
            bool bPointStereo = false;
            if (cosParallaxRays < cosParallaxStereo && cosParallaxRays > 0 && (bStereo1 || bStereo2 || (cosParallaxRays < 0.9996 && mbInertial) || (cosParallaxRays < 0.9998 && !mbInertial))) {
                goodProj = GeometricTools::Triangulate(xn1, xn2, eigTcw1, eigTcw2, x3D);
                if (!goodProj)
                    continue;
            } else if (bStereo1 && cosParallaxStereo1 < cosParallaxStereo2) {
                countStereoAttempt++;
                bPointStereo = true;
                goodProj = mpCurrentKeyFrame->UnprojectStereo(idx1, x3D);
            } else if (bStereo2 && cosParallaxStereo2 < cosParallaxStereo1) {
                countStereoAttempt++;
                bPointStereo = true;
                goodProj = pKF2->UnprojectStereo(idx2, x3D);
            } else {
                continue; // No stereo and very low parallax
            }

            if (goodProj && bPointStereo)
                countStereoGoodProj++;

            if (!goodProj)
                continue;

            // Check triangulation in front of cameras
            float z1 = Rcw1.row(2).dot(x3D) + tcw1(2);
            if (z1 <= 0)
                continue;

            float z2 = Rcw2.row(2).dot(x3D) + tcw2(2);
            if (z2 <= 0)
                continue;

            // Check reprojection error in first keyframe
            const float& sigmaSquare1 = mpCurrentKeyFrame->mvLevelSigma2[kp1.octave];
            const float x1 = Rcw1.row(0).dot(x3D) + tcw1(0);
            const float y1 = Rcw1.row(1).dot(x3D) + tcw1(1);
            const float invz1 = 1.0 / z1;

            if (!bStereo1) {
                cv::Point2f uv1 = pCamera1->project(cv::Point3f(x1, y1, z1));
                float errX1 = uv1.x - kp1.pt.x;
                float errY1 = uv1.y - kp1.pt.y;

                if ((errX1 * errX1 + errY1 * errY1) > 5.991 * sigmaSquare1)
                    continue;

            } else {
                float u1 = fx1 * x1 * invz1 + cx1;
                float u1_r = u1 - mpCurrentKeyFrame->mbf * invz1;
                float v1 = fy1 * y1 * invz1 + cy1;
                float errX1 = u1 - kp1.pt.x;
                float errY1 = v1 - kp1.pt.y;
                float errX1_r = u1_r - kp1_ur;
                if ((errX1 * errX1 + errY1 * errY1 + errX1_r * errX1_r) > 7.8 * sigmaSquare1)
                    continue;
            }

            // Check reprojection error in second keyframe
            const float sigmaSquare2 = pKF2->mvLevelSigma2[kp2.octave];
            const float x2 = Rcw2.row(0).dot(x3D) + tcw2(0);
            const float y2 = Rcw2.row(1).dot(x3D) + tcw2(1);
            const float invz2 = 1.0 / z2;
            if (!bStereo2) {
                cv::Point2f uv2 = pCamera2->project(cv::Point3f(x2, y2, z2));
                float errX2 = uv2.x - kp2.pt.x;
                float errY2 = uv2.y - kp2.pt.y;
                if ((errX2 * errX2 + errY2 * errY2) > 5.991 * sigmaSquare2)
                    continue;
            } else {
                float u2 = fx2 * x2 * invz2 + cx2;
                float u2_r = u2 - mpCurrentKeyFrame->mbf * invz2;
                float v2 = fy2 * y2 * invz2 + cy2;
                float errX2 = u2 - kp2.pt.x;
                float errY2 = v2 - kp2.pt.y;
                float errX2_r = u2_r - kp2_ur;
                if ((errX2 * errX2 + errY2 * errY2 + errX2_r * errX2_r) > 7.8 * sigmaSquare2)
                    continue;
            }

            // Check scale consistency
            Eigen::Vector3f normal1 = x3D - Ow1;
            float dist1 = normal1.norm();

            Eigen::Vector3f normal2 = x3D - Ow2;
            float dist2 = normal2.norm();

            if (dist1 == 0 || dist2 == 0)
                continue;

            if (mbFarPoints && (dist1 >= mThFarPoints || dist2 >= mThFarPoints)) // MODIFICATION
                continue;

            const float ratioDist = dist2 / dist1;
            const float ratioOctave = mpCurrentKeyFrame->mvScaleFactors[kp1.octave] / pKF2->mvScaleFactors[kp2.octave];

            if (ratioDist * ratioFactor < ratioOctave || ratioDist > ratioOctave * ratioFactor)
                continue;

            // Triangulation is succesfull
            MapPoint* pMP = new MapPoint(x3D, mpCurrentKeyFrame, mpAtlas->GetCurrentMap());
            if (bPointStereo)
                countStereo++;

            pMP->AddObservation(mpCurrentKeyFrame, idx1);
            pMP->AddObservation(pKF2, idx2);

            mpCurrentKeyFrame->AddMapPoint(pMP, idx1);
            pKF2->AddMapPoint(pMP, idx2);

            pMP->ComputeDistinctiveDescriptors();

            pMP->UpdateNormalAndDepth();

            mpAtlas->AddMapPoint(pMP);
            mlpRecentAddedMapPoints.push_back(pMP);
        }
    }
}

void LocalMapping::SearchInNeighbors()
{
    // Retrieve neighbor keyframes
    int nn = 10;
    if (mbMonocular)
        nn = 30;
    const vector<KeyFrame*> vpNeighKFs = mpCurrentKeyFrame->GetBestCovisibilityKeyFrames(nn);
    vector<KeyFrame*> vpTargetKFs;
    for (vector<KeyFrame*>::const_iterator vit = vpNeighKFs.begin(), vend = vpNeighKFs.end(); vit != vend; vit++) {
        KeyFrame* pKFi = *vit;
        if (pKFi->isBad() || pKFi->mnFuseTargetForKF == mpCurrentKeyFrame->mnId)
            continue;
        vpTargetKFs.push_back(pKFi);
        pKFi->mnFuseTargetForKF = mpCurrentKeyFrame->mnId;
    }

    // Add some covisible of covisible
    // Extend to some second neighbors if abort is not requested
    for (int i = 0, imax = vpTargetKFs.size(); i < imax; i++) {
        const vector<KeyFrame*> vpSecondNeighKFs = vpTargetKFs[i]->GetBestCovisibilityKeyFrames(20);
        for (vector<KeyFrame*>::const_iterator vit2 = vpSecondNeighKFs.begin(), vend2 = vpSecondNeighKFs.end(); vit2 != vend2; vit2++) {
            KeyFrame* pKFi2 = *vit2;
            if (pKFi2->isBad() || pKFi2->mnFuseTargetForKF == mpCurrentKeyFrame->mnId || pKFi2->mnId == mpCurrentKeyFrame->mnId)
                continue;
            vpTargetKFs.push_back(pKFi2);
            pKFi2->mnFuseTargetForKF = mpCurrentKeyFrame->mnId;
        }
        if (mbAbortBA)
            break;
    }

    // Extend to temporal neighbors
    if (mbInertial) {
        KeyFrame* pKFi = mpCurrentKeyFrame->mPrevKF;
        while (vpTargetKFs.size() < 20 && pKFi) {
            if (pKFi->isBad() || pKFi->mnFuseTargetForKF == mpCurrentKeyFrame->mnId) {
                pKFi = pKFi->mPrevKF;
                continue;
            }
            vpTargetKFs.push_back(pKFi);
            pKFi->mnFuseTargetForKF = mpCurrentKeyFrame->mnId;
            pKFi = pKFi->mPrevKF;
        }
    }

    // Search matches by projection from current KF in target KFs
    ORBmatcher matcher;
    vector<MapPoint*> vpMapPointMatches = mpCurrentKeyFrame->GetMapPointMatches();
    for (vector<KeyFrame*>::iterator vit = vpTargetKFs.begin(), vend = vpTargetKFs.end(); vit != vend; vit++) {
        KeyFrame* pKFi = *vit;

        matcher.Fuse(pKFi, vpMapPointMatches);
        if (pKFi->NLeft != -1)
            matcher.Fuse(pKFi, vpMapPointMatches, true);
    }

    if (mbAbortBA)
        return;

    // Search matches by projection from target KFs in current KF
    vector<MapPoint*> vpFuseCandidates;
    vpFuseCandidates.reserve(vpTargetKFs.size() * vpMapPointMatches.size());

    for (vector<KeyFrame*>::iterator vitKF = vpTargetKFs.begin(), vendKF = vpTargetKFs.end(); vitKF != vendKF; vitKF++) {
        KeyFrame* pKFi = *vitKF;

        vector<MapPoint*> vpMapPointsKFi = pKFi->GetMapPointMatches();

        for (vector<MapPoint*>::iterator vitMP = vpMapPointsKFi.begin(), vendMP = vpMapPointsKFi.end(); vitMP != vendMP; vitMP++) {
            MapPoint* pMP = *vitMP;
            if (!pMP)
                continue;
            if (pMP->isBad() || pMP->mnFuseCandidateForKF == mpCurrentKeyFrame->mnId)
                continue;
            pMP->mnFuseCandidateForKF = mpCurrentKeyFrame->mnId;
            vpFuseCandidates.push_back(pMP);
        }
    }

    matcher.Fuse(mpCurrentKeyFrame, vpFuseCandidates);
    if (mpCurrentKeyFrame->NLeft != -1)
        matcher.Fuse(mpCurrentKeyFrame, vpFuseCandidates, true);

    // Update points
    vpMapPointMatches = mpCurrentKeyFrame->GetMapPointMatches();
    for (size_t i = 0, iend = vpMapPointMatches.size(); i < iend; i++) {
        MapPoint* pMP = vpMapPointMatches[i];
        if (pMP) {
            if (!pMP->isBad()) {
                pMP->ComputeDistinctiveDescriptors();
                pMP->UpdateNormalAndDepth();
            }
        }
    }

    // Update connections in covisibility graph
    mpCurrentKeyFrame->UpdateConnections();
}

void LocalMapping::RequestStop()
{
    unique_lock<mutex> lock(mMutexStop);
    mbStopRequested = true;
    unique_lock<mutex> lock2(mMutexNewKFs);
    mbAbortBA = true;
}

bool LocalMapping::Stop()
{
    unique_lock<mutex> lock(mMutexStop);
    if (mbStopRequested && !mbNotStop) {
        mbStopped = true;
        cout << "Local Mapping STOP" << endl;
        return true;
    }

    return false;
}

bool LocalMapping::isStopped()
{
    unique_lock<mutex> lock(mMutexStop);
    return mbStopped;
}

bool LocalMapping::stopRequested()
{
    unique_lock<mutex> lock(mMutexStop);
    return mbStopRequested;
}

void LocalMapping::Release()
{
    unique_lock<mutex> lock(mMutexStop);
    unique_lock<mutex> lock2(mMutexFinish);
    if (mbFinished)
        return;
    mbStopped = false;
    mbStopRequested = false;
    for (list<KeyFrame*>::iterator lit = mlNewKeyFrames.begin(), lend = mlNewKeyFrames.end(); lit != lend; lit++)
        delete *lit;
    mlNewKeyFrames.clear();

    cout << "Local Mapping RELEASE" << endl;
}

bool LocalMapping::AcceptKeyFrames()
{
    unique_lock<mutex> lock(mMutexAccept);
    return mbAcceptKeyFrames;
}

/**
 * @brief 设置"允许接受关键帧"的状态标志
 */
void LocalMapping::SetAcceptKeyFrames(bool flag)
{
    unique_lock<mutex> lock(mMutexAccept);
    mbAcceptKeyFrames = flag;
}

bool LocalMapping::SetNotStop(bool flag)
{
    unique_lock<mutex> lock(mMutexStop);

    if (flag && mbStopped)
        return false;

    mbNotStop = flag;

    return true;
}

void LocalMapping::InterruptBA()
{
    mbAbortBA = true;
}

void LocalMapping::KeyFrameCulling()
{
    // Check redundant keyframes (only local keyframes)
    // A keyframe is considered redundant if the 90% of the MapPoints it sees, are seen
    // in at least other 3 keyframes (in the same or finer scale)
    // We only consider close stereo points
    const int Nd = 21;
    mpCurrentKeyFrame->UpdateBestCovisibles();
    vector<KeyFrame*> vpLocalKeyFrames = mpCurrentKeyFrame->GetVectorCovisibleKeyFrames();

    float redundant_th;
    if (!mbInertial)
        redundant_th = 0.9;
    else if (mbMonocular)
        redundant_th = 0.9;
    else
        redundant_th = 0.5;

    const bool bInitImu = mpAtlas->isImuInitialized();
    int count = 0;

    // Compoute last KF from optimizable window:
    unsigned int last_ID;
    if (mbInertial) {
        int count = 0;
        KeyFrame* aux_KF = mpCurrentKeyFrame;
        while (count < Nd && aux_KF->mPrevKF) {
            aux_KF = aux_KF->mPrevKF;
            count++;
        }
        last_ID = aux_KF->mnId;
    }

    for (vector<KeyFrame*>::iterator vit = vpLocalKeyFrames.begin(), vend = vpLocalKeyFrames.end(); vit != vend; vit++) {
        count++;
        KeyFrame* pKF = *vit;

        if ((pKF->mnId == pKF->GetMap()->GetInitKFid()) || pKF->isBad())
            continue;
        const vector<MapPoint*> vpMapPoints = pKF->GetMapPointMatches();

        int nObs = 3;
        const int thObs = nObs;
        int nRedundantObservations = 0;
        int nMPs = 0;
        for (size_t i = 0, iend = vpMapPoints.size(); i < iend; i++) {
            MapPoint* pMP = vpMapPoints[i];
            if (pMP) {
                if (!pMP->isBad()) {
                    if (!mbMonocular) {
                        if (pKF->mvDepth[i] > pKF->mThDepth || pKF->mvDepth[i] < 0)
                            continue;
                    }

                    nMPs++;
                    if (pMP->Observations() > thObs) {
                        const int& scaleLevel = (pKF->NLeft == -1) ? pKF->mvKeysUn[i].octave
                            : (i < pKF->NLeft)                     ? pKF->mvKeys[i].octave
                                                                   : pKF->mvKeysRight[i].octave;
                        const map<KeyFrame*, tuple<int, int>> observations = pMP->GetObservations();
                        int nObs = 0;
                        for (map<KeyFrame*, tuple<int, int>>::const_iterator mit = observations.begin(), mend = observations.end(); mit != mend; mit++) {
                            KeyFrame* pKFi = mit->first;
                            if (pKFi == pKF)
                                continue;
                            tuple<int, int> indexes = mit->second;
                            int leftIndex = get<0>(indexes), rightIndex = get<1>(indexes);
                            int scaleLeveli = -1;
                            if (pKFi->NLeft == -1)
                                scaleLeveli = pKFi->mvKeysUn[leftIndex].octave;
                            else {
                                if (leftIndex != -1) {
                                    scaleLeveli = pKFi->mvKeys[leftIndex].octave;
                                }
                                if (rightIndex != -1) {
                                    int rightLevel = pKFi->mvKeysRight[rightIndex - pKFi->NLeft].octave;
                                    scaleLeveli = (scaleLeveli == -1 || scaleLeveli > rightLevel) ? rightLevel
                                                                                                  : scaleLeveli;
                                }
                            }

                            if (scaleLeveli <= scaleLevel + 1) {
                                nObs++;
                                if (nObs > thObs)
                                    break;
                            }
                        }
                        if (nObs > thObs) {
                            nRedundantObservations++;
                        }
                    }
                }
            }
        }

        if (nRedundantObservations > redundant_th * nMPs) {
            if (mbInertial) {
                if (mpAtlas->KeyFramesInMap() <= Nd)
                    continue;

                if (pKF->mnId > (mpCurrentKeyFrame->mnId - 2))
                    continue;

                if (pKF->mPrevKF && pKF->mNextKF) {
                    const float t = pKF->mNextKF->mTimeStamp - pKF->mPrevKF->mTimeStamp;

                    if ((bInitImu && (pKF->mnId < last_ID) && t < 3.) || (t < 0.5)) {
                        pKF->mNextKF->mpImuPreintegrated->MergePrevious(pKF->mpImuPreintegrated);
                        pKF->mNextKF->mPrevKF = pKF->mPrevKF;
                        pKF->mPrevKF->mNextKF = pKF->mNextKF;
                        pKF->mNextKF = NULL;
                        pKF->mPrevKF = NULL;
                        pKF->SetBadFlag();
                    } else if (!mpCurrentKeyFrame->GetMap()->GetIniertialBA2() && ((pKF->GetImuPosition() - pKF->mPrevKF->GetImuPosition()).norm() < 0.02) && (t < 3)) {
                        pKF->mNextKF->mpImuPreintegrated->MergePrevious(pKF->mpImuPreintegrated);
                        pKF->mNextKF->mPrevKF = pKF->mPrevKF;
                        pKF->mPrevKF->mNextKF = pKF->mNextKF;
                        pKF->mNextKF = NULL;
                        pKF->mPrevKF = NULL;
                        pKF->SetBadFlag();
                    }
                }
            } else {
                pKF->SetBadFlag();
            }
        }
        if ((count > 20 && mbAbortBA) || count > 100) {
            break;
        }
    }
}

void LocalMapping::RequestReset()
{
    {
        unique_lock<mutex> lock(mMutexReset);
        cout << "LM: Map reset recieved" << endl;
        mbResetRequested = true;
    }
    cout << "LM: Map reset, waiting..." << endl;

    while (1) {
        {
            unique_lock<mutex> lock2(mMutexReset);
            if (!mbResetRequested)
                break;
        }
        usleep(3000);
    }
    cout << "LM: Map reset, Done!!!" << endl;
}

void LocalMapping::RequestResetActiveMap(Map* pMap)
{
    {
        unique_lock<mutex> lock(mMutexReset);
        cout << "LM: Active map reset recieved" << endl;
        mbResetRequestedActiveMap = true;
        mpMapToReset = pMap;
    }
    cout << "LM: Active map reset, waiting..." << endl;

    while (1) {
        {
            unique_lock<mutex> lock2(mMutexReset);
            if (!mbResetRequestedActiveMap)
                break;
        }
        usleep(3000);
    }
    cout << "LM: Active map reset, Done!!!" << endl;
}

void LocalMapping::ResetIfRequested()
{
    bool executed_reset = false;
    {
        unique_lock<mutex> lock(mMutexReset);
        if (mbResetRequested) {
            executed_reset = true;

            cout << "LM: Reseting Atlas in Local Mapping..." << endl;
            mlNewKeyFrames.clear();
            mlpRecentAddedMapPoints.clear();
            mbResetRequested = false;
            mbResetRequestedActiveMap = false;

            // Inertial parameters
            mTinit = 0.f;
            mbNotBA2 = true;
            mbNotBA1 = true;
            mbBadImu = false;

            mIdxInit = 0;

            cout << "LM: End reseting Local Mapping..." << endl;
        }

        if (mbResetRequestedActiveMap) {
            executed_reset = true;
            cout << "LM: Reseting current map in Local Mapping..." << endl;
            mlNewKeyFrames.clear();
            mlpRecentAddedMapPoints.clear();

            // Inertial parameters
            mTinit = 0.f;
            mbNotBA2 = true;
            mbNotBA1 = true;
            mbBadImu = false;

            mbResetRequested = false;
            mbResetRequestedActiveMap = false;
            cout << "LM: End reseting Local Mapping..." << endl;
        }
    }
    if (executed_reset)
        cout << "LM: Reset free the mutex" << endl;
}

void LocalMapping::RequestFinish()
{
    unique_lock<mutex> lock(mMutexFinish);
    mbFinishRequested = true;
}

bool LocalMapping::CheckFinish()
{
    unique_lock<mutex> lock(mMutexFinish);
    return mbFinishRequested;
}

void LocalMapping::SetFinish()
{
    unique_lock<mutex> lock(mMutexFinish);
    mbFinished = true;
    unique_lock<mutex> lock2(mMutexStop);
    mbStopped = true;
}

bool LocalMapping::isFinished()
{
    unique_lock<mutex> lock(mMutexFinish);
    return mbFinished;
}
/*
 *@brief:获得好的imu初始化值：velocities，gravity direction and IMU biases
 *priorG: gyro先验
 *priorA: acc先验
 *bFIBA : 是否进行FullInertialBA
 */
void LocalMapping::InitializeIMU(float priorG, float priorA, bool bFIBA)
{
    if (mbResetRequested)
        return;

    float minTime; // 初始化需要的最小时间间隔
    int nMinKF; // 初始化需要的最少关键帧数
    if (mbMonocular) {
        minTime = 2.0;
        nMinKF = 10;
    } else {
        minTime = 1.0;
        nMinKF = 10;
    }

    if (mpAtlas->KeyFramesInMap() < nMinKF)
        return;

    // Retrieve all keyframe in temporal order
    // Step 1:按时间顺序收集初始化imu使用的KF
    list<KeyFrame*> lpKF;
    KeyFrame* pKF = mpCurrentKeyFrame;
    while (pKF->mPrevKF) {
        lpKF.push_front(pKF);
        pKF = pKF->mPrevKF;
    }
    lpKF.push_front(pKF);
    vector<KeyFrame*> vpKF(lpKF.begin(), lpKF.end());

    if (vpKF.size() < nMinKF)
        return;

    // imu计算初始时间
    mFirstTs = vpKF.front()->mTimeStamp;
    if (mpCurrentKeyFrame->mTimeStamp - mFirstTs < minTime)
        return;

    bInitializing = true;

    while (CheckNewKeyFrames()) {
        ProcessNewKeyFrame();
        vpKF.push_back(mpCurrentKeyFrame);
        lpKF.push_back(mpCurrentKeyFrame);
    }

    const int N = vpKF.size();
    IMU::Bias b(0, 0, 0, 0, 0, 0);

    // Compute and KF velocities mRwg estimation
    // Step 2:估计KF速度和重力方向
    if (!mpCurrentKeyFrame->GetMap()->isImuInitialized()) {
        Eigen::Matrix3f Rwg;
        Eigen::Vector3f dirG;
        dirG.setZero();
        for (vector<KeyFrame*>::iterator itKF = vpKF.begin(); itKF != vpKF.end(); itKF++) {
            if (!(*itKF)->mpImuPreintegrated)
                continue;
            if (!(*itKF)->mPrevKF)
                continue;

            // 预积分中delta_V 用来表示:Rwb_i.transpose()*(V2 - V1 - g*dt),故此处获得 -(V_i - V_0 - (i-0)*(mRwg*gI)*dt)
            // 应该使用将速度偏差在此处忽略或当做噪声，因为后面会优化mRwg
            dirG -= (*itKF)->mPrevKF->GetImuRotation() * (*itKF)->mpImuPreintegrated->GetUpdatedDeltaVelocity();
            Eigen::Vector3f _vel = ((*itKF)->GetImuPosition() - (*itKF)->mPrevKF->GetImuPosition()) / (*itKF)->mpImuPreintegrated->dT;
            (*itKF)->SetVelocity(_vel);
            (*itKF)->mPrevKF->SetVelocity(_vel);
        }

        // Step 2.1:计算重力方向(与z轴偏差)，用轴角方式表示偏差
        dirG = dirG / dirG.norm();
        Eigen::Vector3f gI(0.0f, 0.0f, -1.0f); //沿-z的归一化的重力数值
        // dirG和gI的模长都是1,故cross为sin，dot为cos
        // 计算旋转轴
        Eigen::Vector3f v = gI.cross(dirG);
        const float nv = v.norm();
        // 计算旋转角
        const float cosg = gI.dot(dirG);
        const float ang = acos(cosg);
        // 计算mRwg，与-Z旋转偏差
        Eigen::Vector3f vzg = v * ang / nv;
        Rwg = Sophus::SO3f::exp(vzg).matrix();
        mRwg = Rwg.cast<double>();
        mTinit = mpCurrentKeyFrame->mTimeStamp - mFirstTs;
    } else {
        mRwg = Eigen::Matrix3d::Identity();
        mbg = mpCurrentKeyFrame->GetGyroBias().cast<double>();
        mba = mpCurrentKeyFrame->GetAccBias().cast<double>();
    }

    mScale = 1.0;

    mInitTime = mpTracker->mLastFrame.mTimeStamp - vpKF.front()->mTimeStamp;

    // Step 3:进行惯性优化
    std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
    // 使用camera初始地图frame的pose与预积分的差值优化
    Optimizer::InertialOptimization(mpAtlas->GetCurrentMap(), mRwg, mScale, mbg, mba, mbMonocular, infoInertial, false, false, priorG, priorA);

    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();

    if (mScale < 1e-1) {
        cout << "scale too small" << endl;
        bInitializing = false;
        return;
    }

    // Before this line we are not changing the map
    // 上面的程序没有改变地图，下面会对地图进行修改
    {
        unique_lock<mutex> lock(mpAtlas->GetCurrentMap()->mMutexMapUpdate);
        if ((fabs(mScale - 1.f) > 0.00001) || !mbMonocular) {
            Sophus::SE3f Twg(mRwg.cast<float>().transpose(), Eigen::Vector3f::Zero());
            mpAtlas->GetCurrentMap()->ApplyScaledRotation(Twg, mScale, true);
            mpTracker->UpdateFrameIMU(mScale, vpKF[0]->GetImuBias(), mpCurrentKeyFrame);
        }

        // Check if initialization OK
        // Step 4:更新关键帧中imu状态
        if (!mpAtlas->isImuInitialized())
            for (int i = 0; i < N; i++) {
                KeyFrame* pKF2 = vpKF[i];
                pKF2->bImu = true;
            }
    }

    mpTracker->UpdateFrameIMU(1.0, vpKF[0]->GetImuBias(), mpCurrentKeyFrame);
    if (!mpAtlas->isImuInitialized()) {
        mpAtlas->SetImuInitialized();
        mpTracker->t0IMU = mpTracker->mCurrentFrame.mTimeStamp;
        mpCurrentKeyFrame->bImu = true;
    }

    // Step 5: 进行完全惯性优化(包括MapPoints)
    std::chrono::steady_clock::time_point t4 = std::chrono::steady_clock::now();
    if (bFIBA) {
        if (priorA != 0.f)
            Optimizer::FullInertialBA(mpAtlas->GetCurrentMap(), 100, false, mpCurrentKeyFrame->mnId, NULL, true, priorG, priorA);
        else
            Optimizer::FullInertialBA(mpAtlas->GetCurrentMap(), 100, false, mpCurrentKeyFrame->mnId, NULL, false);
    }

    std::chrono::steady_clock::time_point t5 = std::chrono::steady_clock::now();

    Verbose::PrintMess("Global Bundle Adjustment finished\nUpdating map ...", Verbose::VERBOSITY_NORMAL);

    // Get Map Mutex
    unique_lock<mutex> lock(mpAtlas->GetCurrentMap()->mMutexMapUpdate);

    unsigned long GBAid = mpCurrentKeyFrame->mnId;

    // Process keyframes in the queue
    while (CheckNewKeyFrames()) {
        ProcessNewKeyFrame();
        vpKF.push_back(mpCurrentKeyFrame);
        lpKF.push_back(mpCurrentKeyFrame);
    }

    // Correct keyframes starting at map first keyframe
    list<KeyFrame*> lpKFtoCheck(mpAtlas->GetCurrentMap()->mvpKeyFrameOrigins.begin(), mpAtlas->GetCurrentMap()->mvpKeyFrameOrigins.end());

    while (!lpKFtoCheck.empty()) {
        KeyFrame* pKF = lpKFtoCheck.front();
        const set<KeyFrame*> sChilds = pKF->GetChilds();
        Sophus::SE3f Twc = pKF->GetPoseInverse();
        for (set<KeyFrame*>::const_iterator sit = sChilds.begin(); sit != sChilds.end(); sit++) {
            KeyFrame* pChild = *sit;
            if (!pChild || pChild->isBad())
                continue;

            if (pChild->mnBAGlobalForKF != GBAid) {
                Sophus::SE3f Tchildc = pChild->GetPose() * Twc;
                pChild->mTcwGBA = Tchildc * pKF->mTcwGBA;

                Sophus::SO3f Rcor = pChild->mTcwGBA.so3().inverse() * pChild->GetPose().so3();
                if (pChild->isVelocitySet()) {
                    pChild->mVwbGBA = Rcor * pChild->GetVelocity();
                } else {
                    Verbose::PrintMess("Child velocity empty!! ", Verbose::VERBOSITY_NORMAL);
                }

                pChild->mBiasGBA = pChild->GetImuBias();
                pChild->mnBAGlobalForKF = GBAid;
            }
            lpKFtoCheck.push_back(pChild);
        }

        pKF->mTcwBefGBA = pKF->GetPose();
        pKF->SetPose(pKF->mTcwGBA);

        if (pKF->bImu) {
            pKF->mVwbBefGBA = pKF->GetVelocity();
            pKF->SetVelocity(pKF->mVwbGBA);
            pKF->SetNewBias(pKF->mBiasGBA);
        } else {
            cout << "KF " << pKF->mnId << " not set to inertial!! \n";
        }

        lpKFtoCheck.pop_front();
    }

    // Correct MapPoints
    const vector<MapPoint*> vpMPs = mpAtlas->GetCurrentMap()->GetAllMapPoints();

    for (size_t i = 0; i < vpMPs.size(); i++) {
        MapPoint* pMP = vpMPs[i];

        if (pMP->isBad())
            continue;

        if (pMP->mnBAGlobalForKF == GBAid) {
            // If optimized by Global BA, just update
            pMP->SetWorldPos(pMP->mPosGBA);
        } else {
            // Update according to the correction of its reference keyframe
            KeyFrame* pRefKF = pMP->GetReferenceKeyFrame();

            if (pRefKF->mnBAGlobalForKF != GBAid)
                continue;

            // Map to non-corrected camera
            Eigen::Vector3f Xc = pRefKF->mTcwBefGBA * pMP->GetWorldPos();

            // Backproject using corrected camera
            pMP->SetWorldPos(pRefKF->GetPoseInverse() * Xc);
        }
    }

    Verbose::PrintMess("Map updated!", Verbose::VERBOSITY_NORMAL);

    mnKFs = vpKF.size();
    mIdxInit++;

    for (list<KeyFrame*>::iterator lit = mlNewKeyFrames.begin(), lend = mlNewKeyFrames.end(); lit != lend; lit++) {
        (*lit)->SetBadFlag();
        delete *lit;
    }
    mlNewKeyFrames.clear();

    mpTracker->mState = Tracking::OK;
    bInitializing = false;

    mpCurrentKeyFrame->GetMap()->IncreaseChangeIndex();

    return;
}

void LocalMapping::ScaleRefinement()
{
    // Minimum number of keyframes to compute a solution
    // Minimum time (seconds) between first and last keyframe to compute a solution. Make the difference between monocular and stereo
    // unique_lock<mutex> lock0(mMutexImuInit);
    if (mbResetRequested)
        return;

    // Retrieve all keyframes in temporal order
    list<KeyFrame*> lpKF;
    KeyFrame* pKF = mpCurrentKeyFrame;
    while (pKF->mPrevKF) {
        lpKF.push_front(pKF);
        pKF = pKF->mPrevKF;
    }
    lpKF.push_front(pKF);
    vector<KeyFrame*> vpKF(lpKF.begin(), lpKF.end());

    while (CheckNewKeyFrames()) {
        ProcessNewKeyFrame();
        vpKF.push_back(mpCurrentKeyFrame);
        lpKF.push_back(mpCurrentKeyFrame);
    }

    const int N = vpKF.size();

    mRwg = Eigen::Matrix3d::Identity();
    mScale = 1.0;

    std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
    Optimizer::InertialOptimization(mpAtlas->GetCurrentMap(), mRwg, mScale);
    std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();

    if (mScale < 1e-1) // 1e-1
    {
        cout << "scale too small" << endl;
        bInitializing = false;
        return;
    }

    Sophus::SO3d so3wg(mRwg);
    // Before this line we are not changing the map
    unique_lock<mutex> lock(mpAtlas->GetCurrentMap()->mMutexMapUpdate);
    std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
    if ((fabs(mScale - 1.f) > 0.002) || !mbMonocular) {
        Sophus::SE3f Tgw(mRwg.cast<float>().transpose(), Eigen::Vector3f::Zero());
        mpAtlas->GetCurrentMap()->ApplyScaledRotation(Tgw, mScale, true);
        mpTracker->UpdateFrameIMU(mScale, mpCurrentKeyFrame->GetImuBias(), mpCurrentKeyFrame);
    }
    std::chrono::steady_clock::time_point t3 = std::chrono::steady_clock::now();

    for (list<KeyFrame*>::iterator lit = mlNewKeyFrames.begin(), lend = mlNewKeyFrames.end(); lit != lend; lit++) {
        (*lit)->SetBadFlag();
        delete *lit;
    }
    mlNewKeyFrames.clear();

    double t_inertial_only = std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0).count();

    // To perform pose-inertial opt w.r.t. last keyframe
    mpCurrentKeyFrame->GetMap()->IncreaseChangeIndex();

    return;
}

bool LocalMapping::IsInitializing()
{
    return bInitializing;
}

double LocalMapping::GetCurrKFTime()
{

    if (mpCurrentKeyFrame) {
        return mpCurrentKeyFrame->mTimeStamp;
    } else
        return 0.0;
}

KeyFrame* LocalMapping::GetCurrKF()
{
    return mpCurrentKeyFrame;
}

} // namespace ORB_SLAM
