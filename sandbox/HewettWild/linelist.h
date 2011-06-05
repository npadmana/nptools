#ifndef LINELIST_H
#define LINELIST_H

#include <map>
#include <string>
#include <utility>
#include <Eigen/Dense>

class LineList
{
private :
    // Store the line list in a map
    typedef std::multimap<std::string, double> linemap;
    linemap llist;
public:
    LineList();
    int nlines(std::string linestr);
    Eigen::VectorXd wave(std::string linestr);
};

#endif // LINELIST_H
