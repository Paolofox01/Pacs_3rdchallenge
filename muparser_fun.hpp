#ifndef MUPARSERFUN_H
#define MUPARSERFUN_H


#include <muParser.h>
#include <Eigen/Eigen>
#include <memory>
#include <string>
#include <cmath>

class MuparserFun
{
public:
  MuparserFun(const MuparserFun &m)
    : m_parser(m.m_parser)
  {
    m_var.resize(2);
    m_parser.DefineVar("x", &m_var[0]);
    m_parser.DefineVar("y", &m_var[1]);
  };

  MuparserFun(const std::string &s)
  {
    try
      {
        m_var.resize(2);
        m_parser.DefineVar("x", &m_var[0]);
        m_parser.DefineVar("y", &m_var[1]);
        m_parser.SetExpr(s);
      }
    catch (mu::Parser::exception_type &e)
      {
        std::cerr << e.GetMsg() << std::endl;
      }
  };

  double
  operator()(const double &x, const double &y)
  {
    m_var << x, y;
    double z = m_parser.Eval();
    
    return z;
  };

  Eigen::VectorXd GetVar()
  {
      return m_var;
  }

  
  mu::Parser Getparser()
  {
      return m_parser;
  }


private:
  Eigen::VectorXd m_var;
  mu::Parser m_parser;
};


#endif