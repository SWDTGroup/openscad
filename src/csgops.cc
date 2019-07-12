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

#include "csgnode.h"

#include "evalcontext.h"
#include "module.h"
#include "csgterm.h"
#include "builtin.h"
#include <sstream>
#include <assert.h>

#include "evalcontext.h"
#include <vector>
#include <boost/assign/std/vector.hpp>
using namespace boost::assign; //add by Look bring 'operator+=()' into scope

class CsgModule : public AbstractModule
{
public:
	OpenSCADOperator type;
	CsgModule(OpenSCADOperator type) : type(type) { }
	virtual AbstractNode *instantiate(const Context *ctx, const ModuleInstantiation *inst, EvalContext *evalctx) const;
};

AbstractNode *CsgModule::instantiate(const Context* ctx, const ModuleInstantiation *inst, EvalContext *evalctx) const
{
	inst->scope.apply(*evalctx);
	CsgNode *node = new CsgNode(inst, type);
	//add by Look begin
	if (type == OPENSCAD_CARVE)
	{
		AssignmentList args;
		args += Assignment("depth");
		args += Assignment("subdivideLen");
		Context c(ctx);
		c.setVariables(args, evalctx);
		inst->scope.apply(*evalctx);

		ValuePtr depth = c.lookup_variable("depth");
		if (depth->type() == Value::NUMBER)
		{
			node->carve_depth  = depth->toDouble();
		}
		ValuePtr subdivideLen = c.lookup_variable("subdivideLen");
		if (subdivideLen->type() == Value::NUMBER)
		{
			node->carve_subdivideLen  = subdivideLen->toDouble();
		}

		//node->fn = c.lookup_variable("$fn")->toDouble();
		//node->fs = c.lookup_variable("$fs")->toDouble();
		//node->fa = c.lookup_variable("$fa")->toDouble();
	}
	//add by Look end
	std::vector<AbstractNode *> instantiatednodes = inst->instantiateChildren(evalctx);
	node->children.insert(node->children.end(), instantiatednodes.begin(), instantiatednodes.end());
	return node;
}

std::string CsgNode::toString() const
{
	std::stringstream stream;
	stream << this->name() << "(depth=" << carve_depth << ", subdivideLen=" << carve_subdivideLen << ")";
	return stream.str();
}

std::string CsgNode::name() const
{
	switch (this->type) {
	case OPENSCAD_CARVE:
		return "carve";
		break;
	case OPENSCAD_UNION:
		return "union";
		break;
	case OPENSCAD_DIFFERENCE:
		return "difference";
		break;
	case OPENSCAD_INTERSECTION:
		return "intersection";
		break;
	case OPENSCAD_APPEND:
		return "append";
		break;
	default:
		assert(false);
	}
	return "internal_error";
}

void register_builtin_csgops()
{
	Builtins::init("carve", new CsgModule(OPENSCAD_CARVE));	//add by Look
	Builtins::init("union", new CsgModule(OPENSCAD_UNION));
	Builtins::init("difference", new CsgModule(OPENSCAD_DIFFERENCE));
	Builtins::init("intersection", new CsgModule(OPENSCAD_INTERSECTION));
	Builtins::init("append", new CsgModule(OPENSCAD_APPEND));
}

