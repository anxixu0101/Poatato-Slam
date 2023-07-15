
/**
 * @file icp.h
 * @author virtual虚函数
 * @brief   头文件
 * @version 0.1
 * @date 2022-04-27
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include "Eigen/Eigen"
#include <vector>

#ifndef ICP_H
#define ICP_H


class Icp
{

public:

typedef struct{
    Eigen::Matrix4d trans;
    std::vector<float> distances;
    
} IcpOut; //Icp的输出结果

typedef struct{
    std::vector<float> distances;
    std::vector<int> indices;
} NeighBor;

Eigen::Matrix4d bestFitTransform(const Eigen::MatrixXd &source, const Eigen::MatrixXd &target);

IcpOut icp(const Eigen::MatrixXd &source, const Eigen::MatrixXd &target, int max_iterations, int tolerance );


NeighBor nearestNeighbot(const Eigen::MatrixXd &src, const Eigen::MatrixXd &dst);
float getDistance(const Eigen::Vector3d &source, const Eigen::Vector3d &target);

private:

};
#endif