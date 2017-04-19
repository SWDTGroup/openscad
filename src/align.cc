#include "alignNode.h"

#include "module.h"
#include "evalcontext.h"
#include "printutils.h"
#include "builtin.h"
#include "visitor.h"
#include "polyset.h"

#include <assert.h>
#include <sstream>
#include <boost/assign/std/vector.hpp>
using namespace boost::assign; // bring 'operator+=()' into scope

class AlignModule : public AbstractModule
{
public:
	AlignModule () { }
	virtual AbstractNode *instantiate(const Context *ctx, const ModuleInstantiation *inst, EvalContext *evalctx) const;
};

AbstractNode *AlignModule::instantiate(const Context *ctx, const ModuleInstantiation *inst, EvalContext *evalctx) const
{
	AlignNode*node = new AlignNode(inst);

	AssignmentList args;
	args += Assignment("plane");

	Context c(ctx);
	c.setVariables(args, evalctx);
	inst->scope.apply(*evalctx);


	ValuePtr ns = c.lookup_variable("plane");

	if ( ns->type() == Value::VECTOR) {
		const Value::VectorType &vs = ns->toVector();
		if ( vs.size() >= 3 ) 
		{
			node->m_plane[0] = vs[0]->toDouble();
			node->m_plane[1] = vs[1]->toDouble();
			node->m_plane[2] = vs[2]->toDouble();
		}
	}


	std::vector<AbstractNode *> instantiatednodes = inst->instantiateChildren(evalctx);
	node->children.insert(node->children.end(), instantiatednodes.begin(), instantiatednodes.end());

	return node;
}

std::string AlignNode::toString() const
{
	std::stringstream stream;

	stream << "Align(plane = ["<< m_plane[0] << "," <<  m_plane[1] <<"," << m_plane[2]<< "])";

	return stream.str();
}

void register_builtin_align()
{
	Builtins::init("align", new AlignModule());
}
