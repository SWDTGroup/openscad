//add by Look
#pragma once

#include "value.h"
#include "node.h"
#include "visitor.h"
#include "linalg.h"

class PolarizationNode : public AbstractNode
{
public:
	PolarizationNode(const ModuleInstantiation *mi) : AbstractNode(mi) { 
	fs = 0;
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
	double fs; //>0��ζ�з֣�����x���з֣����len=fs
	Value::VectorType params;
	
	std::string exp; 

};
