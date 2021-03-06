/*
 *  OpenSCAD (www.openscad.org)
 *  Copyright (C) 2009-2011 Clifford Wolf <clifford@clifford.at> and
 *                          Marius Kintel <marius@kintel.net>
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  As a special exception, you have permission to link this program
 *  with the CGAL library and distribute executables, as long as you
 *  follow the requirements of the GNU GPL in regard to all of the
 *  software in the executable aside from CGAL.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */
 //add by Look

#include "polarizationNode.h"
#include "module.h"
#include "evalcontext.h"
#include "polyset.h"
#include "builtin.h"
#include "value.h"
#include "printutils.h"
#include <sstream>
#include <vector>
#include <assert.h>
#include <boost/assign/std/vector.hpp>
using namespace boost::assign; // bring 'operator+=()' into scope

enum polarization_type_e {
	polarization_normal
};

class PolarizationModule : public AbstractModule
{
public:
	polarization_type_e type;
	PolarizationModule(polarization_type_e type) : type(type) { }
	virtual AbstractNode *instantiate(const Context *ctx, const ModuleInstantiation *inst, EvalContext *evalctx) const;
};

AbstractNode *PolarizationModule::instantiate(const Context *ctx, const ModuleInstantiation *inst, EvalContext *evalctx) const
{
	PolarizationNode *node = new PolarizationNode(inst);

	
	AssignmentList args;

	switch (this->type) {
	case polarization_normal:
		args += Assignment("arr");
		break;
	default:
		assert(false);
	}

	Context c(ctx);
	c.setVariables(args, evalctx);
	inst->scope.apply(*evalctx);

	if (this->type == polarization_normal)
	{
		ValuePtr arr = c.lookup_variable("arr");
		if ( arr->toVector().size() != 4 ) {
			assert(!"polarization's arguments count is not 4!");
		}
		arr->toVector()[0]->getDouble(node->o_size[0]);
		arr->toVector()[1]->getDouble(node->o_size[1]);
		arr->toVector()[2]->getDouble(node->k_xy[0]);
		arr->toVector()[3]->getDouble(node->k_xy[1]);
		// delete output
		//printf("compiled polarization(%f, %f, %f, %f)\n", node->o_size[0], node->o_size[1], node->k_xy[0], node->k_xy[1]);
	}

	std::vector<AbstractNode *> instantiatednodes = inst->instantiateChildren(evalctx);
	node->children.insert(node->children.end(), instantiatednodes.begin(), instantiatednodes.end());

	return node;
}

std::string PolarizationNode::toString() const
{
	std::stringstream stream;

	stream << "polarization([";
	stream << o_size[0] << ", ";
	stream << o_size[1] << ", ";
	stream << k_xy[0] << ", ";
	stream << k_xy[1];
	stream << "])";

	return stream.str();
}

std::string PolarizationNode::name() const
{
	return "polarization";
}

void register_builtin_polarization()
{
	Builtins::init("polarization", new PolarizationModule(polarization_normal));
}
