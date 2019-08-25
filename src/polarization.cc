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
#include "fileutils.h"

#include <iostream>
#include <fstream>
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
	polarization_normal,
	lua_normal

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
		args += Assignment("angle");
		args += Assignment("max_len");
		break;
	case lua_normal:
		args += Assignment("exp");
		args += Assignment("max_len");
		args += Assignment("params");
		break;
	default:
		assert(false);
	}

	Context c(ctx);
	c.setVariables(args, evalctx);
	inst->scope.apply(*evalctx);

	if (this->type == polarization_normal)
	{
		ValuePtr angle = c.lookup_variable("angle");
		if (angle->type() == Value::NUMBER)
		{
			node->angle  = angle->toDouble();
		}
	}
	else
	{
		ValuePtr v = c.lookup_variable("file");

		std::string filename = lookup_file(v->isUndefined() ? "" : v->toString(), inst->path(), ctx->documentPath());

		   std::ifstream inFile(filename.c_str(), std::ios::in | std::ios::binary);
		    std::ostringstream oss;
		    oss << inFile.rdbuf();
		    node->exp = oss.str();
		    inFile.close();
		    
		ValuePtr params = c.lookup_variable("params");
		if (!params->isUndefined()) 
			node->params = params->toVector();

	}
	node->fs = c.lookup_variable("max_len")->toDouble();


	std::vector<AbstractNode *> instantiatednodes = inst->instantiateChildren(evalctx);
	node->children.insert(node->children.end(), instantiatednodes.begin(), instantiatednodes.end());

	return node;
}

std::string PolarizationNode::toString() const
{
	std::stringstream stream;

	stream << "polarization(angle = " << angle << ", " ;
	stream  <<  "max_len = " << this->fs << ")";

	return stream.str();
}

std::string PolarizationNode::name() const
{
	return "polarization";
}

void register_builtin_polarization()
{
	Builtins::init("sz_polarization", new PolarizationModule(polarization_normal));
	Builtins::init("sz_lua", new PolarizationModule(lua_normal));
}
