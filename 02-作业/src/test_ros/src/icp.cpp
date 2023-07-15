/**
 * @file icp.cpp
 * @author virtual
 * @brief
 * @version 0.1
 * @date 2022-04-26
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <iostream>
#include <numeric>
#include "../include/icp.hpp"
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

/**
 * @brief  ICP算法利用SVD得到R和T的关键算法步骤
 * 
 * @param source 
 * @param target :source 和 target的数量是一定相等的，因为target是根据source最近点选中的点云
 * @return Eigen::Matrix4d 
 */

Eigen::Matrix4d Icp::bestFitTransform(const Eigen::MatrixXd &source, const Eigen::MatrixXd &target)
{

    Eigen::Matrix4d T = Eigen::MatrixXd::Identity(4, 4);
    Eigen::Vector3d centroid_source(0, 0, 0); // 点云中心
    Eigen::Vector3d centroid_target(0, 0, 0);
    Eigen::MatrixXd source_point = source;
    Eigen::MatrixXd target_point = target;
    int row_source = source.rows(); // 得到原始点云的行数（只有3列）
    int row_target = target.rows();

    for (int i = 0; i < row_source; i++)
    {
        centroid_source += source.block<1, 3>(i, 0).transpose();
    }
    for (int i = 0; i < row_target; i++)
    {
        centroid_target += target.block<1, 3>(i, 0).transpose();
    }

    centroid_source /= row_source;
    centroid_target /= row_target; // 分别求出点云中心

    // 去中心化后的点云
    for (int i = 0; i < row_source; i++)
    {
        source_point.block<1, 3>(i, 0) = source.block<1, 3>(i, 0) - centroid_source.transpose();
    }
    for (int i = 0; i < row_target; i++)
    {
        target_point.block<1, 3>(i, 0) = target.block<1, 3>(i, 0) - centroid_target.transpose();
    }

    Eigen::MatrixXd H = source_point.transpose() * target_point;
    Eigen::MatrixXd U;
    Eigen::VectorXd S;
    Eigen::MatrixXd V;
    Eigen::MatrixXd Vt;
    Eigen::Matrix3d R;
    Eigen::Vector3d t;

    JacobiSVD<Eigen::MatrixXd> svd(H, ComputeFullU | ComputeFullV);
    U = svd.matrixU();
    S = svd.singularValues();
    V = svd.matrixV();
    Vt = V.transpose();

    R = Vt.transpose() * U.transpose();

    if (R.determinant() < 0)
    {
        Vt.block<1, 3>(2, 0) *= -1;
        R = Vt.transpose() * U.transpose();
    }

    t = centroid_target - R * centroid_source;

    T.block<3, 3>(0, 0) = R;
    T.block<3, 1>(0, 3) = t;
    return T;
}

/**
 * @brief 
 * 
 * @param A 
 * @param B 
 * @param max_iterations 
 * @param tolerance 
 * @return Icp::IcpOut 
 */

Icp::IcpOut Icp::icp(const Eigen::MatrixXd &source, const Eigen::MatrixXd &target, int max_iterations, int tolerance){
    int row = source.rows();
    Eigen::MatrixXd src = Eigen::MatrixXd::Ones(3+1,row);
    Eigen::MatrixXd src3d = Eigen::MatrixXd::Ones(3,row);
    Eigen::MatrixXd dst = Eigen::MatrixXd::Ones(3+1,row);
    NeighBor neighbor; 
    Eigen::Matrix4d T;
    Eigen::MatrixXd dst_chorder = Eigen::MatrixXd::Ones(3,row);
    IcpOut result;
    int iter = 0;

    for (int i = 0; i<row; i++){
        src.block<3,1>(0,i) = source.block<1,3>(i,0).transpose();
        src3d.block<3,1>(0,i) = source.block<1,3>(i,0).transpose();
        dst.block<3,1>(0,i) = target.block<1,3>(i,0).transpose();

    }

    double prev_error = 0;
    double mean_error = 0;
    for (int i=0; i<max_iterations; i++){
        neighbor = nearestNeighbot(src3d.transpose(),target);

        for(int j=0; j<row; j++){
            dst_chorder.block<3,1>(0,j) = dst.block<3,1>(0,neighbor.indices[j]);
        }

        T = bestFitTransform(src3d.transpose(),dst_chorder.transpose());

        src = T*src;
        for(int j=0; j<row; j++){
            src3d.block<3,1>(0,j) = src.block<3,1>(0,j);
        }

        mean_error = std::accumulate(neighbor.distances.begin(),neighbor.distances.end(),0.0)/neighbor.distances.size();
        if (abs(prev_error - mean_error) < tolerance){
            break;
        }
        prev_error = mean_error;
        iter = i+2;
    }

    T = bestFitTransform(source,src3d.transpose());
    result.trans = T;
    result.distances = neighbor.distances;
  

    return result;
}



/**
 * @brief  得到基于原始点云的最近点云对
 * 
 * @param src 原始点云
 * @param target 目标点云
 * @return NeighBor 
 */
Icp::NeighBor Icp::nearestNeighbot(const Eigen::MatrixXd &src, const Eigen::MatrixXd &target){
    int row_src = src.rows();
    int row_dst = target.rows();
    Eigen::Vector3d vec_src;
    Eigen::Vector3d vec_target;
    NeighBor neigh_point;
    float min = 100;
    int index = 0;
    float dist_temp = 0;

    for(int i=0; i < row_src; i++){
        vec_src = src.block<1,3>(i,0).transpose();
        min = 100;
        index = 0;
        dist_temp = 0;
        for(int j=0; j < row_dst; j++){
            vec_target = target.block<1,3>(j,0).transpose();
            dist_temp = getDistance(vec_src,vec_target);
            if (dist_temp < min){
                min = dist_temp;
                index = j;
            }
        }
        // cout << min << " " << index << endl;
        // neigh.distances[ii] = min;
        // neigh.indices[ii] = index;
        neigh_point.distances.push_back(min);
        neigh_point.indices.push_back(index);
    }

    return neigh_point;
}

float Icp::getDistance(const Eigen::Vector3d &source, const Eigen::Vector3d &target){
    return sqrt((source[0]-target[0])*(source[0]-target[0]) + (source[1]-target[1])*(source[1]-target[1]) + (source[2]-target[2])*(source[2]-target[2]));
}
