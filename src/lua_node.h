#pragma once

#include "node.h"
#include "visitor.h"
#include "value.h"

class LuaNode : public AbstractPolyNode
{
public:
	LuaNode(const ModuleInstantiation *mi) : AbstractPolyNode(mi) {
	}
  virtual Response accept(class State &state, Visitor &visitor) const {
		return visitor.visit(state, *this);
	}
	virtual std::string toString() const;
	virtual std::string name() const { return "sz_lua2"; }

	Filename filename;

};
