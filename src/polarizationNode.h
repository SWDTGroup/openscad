//add by Look
#pragma once

#include "node.h"
#include "NodeVisitor.h"
#include "linalg.h"

class PolarizationNode : public AbstractNode
{
public:
	VISITABLE();
	PolarizationNode(const ModuleInstantiation *mi) : AbstractNode(mi) { }
  virtual Response accept(class State &state, Visitor &visitor) const {
		return visitor.visit(state, *this);
	}
	virtual std::string toString() const;
	virtual std::string name() const;

	double o_size[2];
	double k_xy[2];
};
