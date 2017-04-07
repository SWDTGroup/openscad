#pragma once

#include "node.h"
#include "visitor.h"
#include "value.h"
#include "clipper-utils.h"

class DecimationNode : public AbstractPolyNode
{
public:
	DecimationNode(const ModuleInstantiation *mi) : AbstractPolyNode(mi), target(0) { }
        virtual Response accept(class State &state, Visitor &visitor) const {
		return visitor.visit(state, *this);
	}
	virtual std::string toString() const;
	virtual std::string name() const { return "decimate"; }

    unsigned int target;
};
