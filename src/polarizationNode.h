//add by Look
#pragma once

#include "node.h"
#include "visitor.h"
#include "linalg.h"

class PolarizationNode : public AbstractNode
{
public:
	PolarizationNode(const ModuleInstantiation *mi) : AbstractNode(mi) { 
	fn = fs = fa = 0;
	o_size = 0;
	angle = 360;
}
  virtual Response accept(class State &state, Visitor &visitor) const {
		return visitor.visit(state, *this);
	}
	virtual std::string toString() const;
	virtual std::string name() const;

	double o_size;
	double angle;
	double fn, fs, fa;
	std::string exp; 

};
