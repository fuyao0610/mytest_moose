/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/
#ifndef MMSCONVECTION_H_
#define MMSCONVECTION_H_

#include "Kernel.h"

class MMSConvection;

template<>
InputParameters validParams<MMSConvection>();

class MMSConvection : public Kernel
{
public:
  MMSConvection(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

  RealVectorValue velocity;

  Real _x;
  Real _y;
  Real _z;

};

#endif //MMSCONVECTION_H_
