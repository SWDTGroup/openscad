#pragma once

#include "node.h"
#include "linalg.h"
#include "visitor.h"
#include <string>

class AlignNode : public AbstractNode
{
public:
	AlignNode(const ModuleInstantiation *mi) : AbstractNode(mi) {
		m_plane = Vector3d(0,0,1);
	}
  virtual Response accept(class State &state, Visitor &visitor) const {
		return visitor.visit(state, *this);
	}
	virtual std::string toString() const;
	virtual std::string name() const { return "align"; }

	Vector3d m_plane;
};
