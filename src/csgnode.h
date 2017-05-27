#pragma once

#include "node.h"
#include "visitor.h"
#include "enums.h"
#include "linalg.h"	//add by Look

class CsgNode : public AbstractNode
{
public:
	OpenSCADOperator type;
	CsgNode(const ModuleInstantiation *mi, OpenSCADOperator type) : AbstractNode(mi), type(type) {
		carve_depth = 0.4;		//add by Look
		carve_subdivideLen = 0.2;	//add by Look
	}
  virtual Response accept(class State &state, Visitor &visitor) const {
		return visitor.visit(state, *this);
	}
	virtual std::string toString() const;
	virtual std::string name() const;
	double carve_depth, carve_subdivideLen;	//for carve() add by Look
};
