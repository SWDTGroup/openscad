
#include "lua_node.h"

#include "module.h"
#include "evalcontext.h"
#include "printutils.h"
#include "fileutils.h"
#include "builtin.h"
#include "calc.h"
#include "polyset.h"
#include "mathc99.h" 

#include <sstream>
#include <boost/assign/std/vector.hpp>
using namespace boost::assign; // bring 'operator+=()' into scope

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

class LuaModule : public AbstractModule
{
public:
	LuaModule() { }
	virtual AbstractNode *instantiate(const Context *ctx, const ModuleInstantiation *inst, EvalContext *evalctx) const;
};

AbstractNode *LuaModule::instantiate(const Context *ctx, const ModuleInstantiation *inst, EvalContext *evalctx) const
{
	LuaNode *node = new LuaNode(inst);

	AssignmentList args;
	args += Assignment("file");

	Context c(ctx);
	c.setVariables(args, evalctx);
	inst->scope.apply(*evalctx);


	ValuePtr file = c.lookup_variable("file");
	

	if (!file->isUndefined() && file->type() == Value::STRING) {
		node->filename = lookup_file(file->toString(), inst->path(), c.documentPath());
	}

	std::vector<AbstractNode *> instantiatednodes = inst->instantiateChildren(evalctx);
	node->children.insert(node->children.end(), instantiatednodes.begin(), instantiatednodes.end());

	return node;
}

std::string LuaNode::toString() const
{
	std::stringstream stream;

	stream << this->name() << "(";
	if (!this->filename.empty()) { // Ignore deprecated parameters if empty 
		fs::path path((std::string)this->filename);
		stream <<
			"file = " << this->filename << ", "
#ifndef OPENSCAD_TESTING
			// timestamp is needed for caching, but disturbs the test framework
			<< "timestamp = " << (fs::exists(path) ? fs::last_write_time(path) : 0) << ", "
#endif
			;
	}

	return stream.str();
}

void register_builtin_sz_lua2()
{

	Builtins::init("sz_lua2", new LuaModule());
}
