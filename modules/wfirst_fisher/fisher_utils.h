/* Basic Fisher matrix manipulations.
 *
 * Nikhil Padmanabhan
 */
#ifndef FISHER_UTILS_H
#define FISHER_UTILS_H

#include <vector>
#include <eigen3/Eigen/Dense>

Eigen::MatrixXd delete_parameter(const Eigen::MatrixXd& fish_in, int par);
Eigen::MatrixXd delete_row(const Eigen::MatrixXd& in, int irow);
Eigen::MatrixXd delete_col(const Eigen::MatrixXd& in, int icol);

Eigen::MatrixXd marginalize(const Eigen::MatrixXd& in, std::vector<int> params);



#endif // FISHER_UTILS_H
