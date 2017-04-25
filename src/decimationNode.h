#pragma once

#include "node.h"
#include "visitor.h"
#include "value.h"
#include "clipper-utils.h"

class DecimationNode : public AbstractPolyNode
{
public:
	enum DecimationOp {DecimationOp_normal, DecimationOp_keepmain};
   	unsigned int target;
	bool keep_main;

	DecimationNode(const ModuleInstantiation *mi) : AbstractPolyNode(mi), target(0) , keep_main(false){ }
        virtual Response accept(class State &state, Visitor &visitor) const {
		return visitor.visit(state, *this);
	}
	virtual std::string toString() const;
	virtual std::string name() const { return "decimate"; }
};
