/*
Copyright 2022 René Ferdinand Rivera Morell
Distributed under the Boost Software License, Version 1.0.
(See accompanying file LICENSE.txt or https://www.bfgroup.xyz/b2/LICENSE.txt)
*/

#include "mod_jam_class.h"

#include "bindjam.h"
#include "compile.h"
#include "modules.h"
#include "variable.h"

#include <unordered_set>

namespace b2 { namespace jam { namespace klass {

std::string make(std::tuple<value_ref, list_ref> name_arg1,
	const lists & rest,
	bind::context_ref_ context_ref)
{
	static std::size_t next_instance = 1;

	std::string id = "object(";
	id += std::get<0>(name_arg1);
	id += ")@" + std::to_string(next_instance);

	std::string type = "class@";
	type += std::get<0>(name_arg1);

	module_ptr instance = bindmodule(value_ref(id));
	module_ptr class_module = bindmodule(value_ref(type));
	instance->class_module = class_module;
	module_set_fixed_variables(instance, class_module->num_fixed_variables);

	variable(id, "__class__") = std::get<0>(name_arg1);
	variable(id, "__name__") = id;

	import_module(*list_ref(id), root_module());

	frame * outer = context_ref.get<jam_context>().frame;
	lists args;
	args | *std::get<1>(name_arg1) | rest;
	list_ref(call_member_rule(
				 value_ref("__init__"), outer, list_ref(id), std::move(args)),
		true);

	// Bump the next unique object name.
	next_instance += 1;

	// Return the name of the new instance.
	return id;
}

list_ref bases(std::string class_name)
{
	return *variable("class@" + class_name, "__bases__");
}

bool is_derived(value_ref class_name, list_cref class_bases)
{
	bool found = false;
	list_ref stack(class_name);
	std::unordered_set<value_ref, value_ref::hash_function,
		value_ref::equal_function>
		visited;
	while (!found && !stack.empty())
	{
		value_ref top = stack[0];
		stack.pop_front();
		if (visited.count(top) == 0)
		{
			visited.insert(top);
			stack.append(bases(top));

			list_cref::size_type bases_found = 0;
			for (auto class_base : class_bases)
				if (visited.count(class_base) == 0)
					break;
				else
					bases_found += 1;
			found = bases_found == class_bases.length();
		}
	}
	return found;
}

bool is_instance(std::string value)
{
	bool has_prefix = value.compare(0, 7, "object(") == 0;
	bool has_mid = value.find(")@", 7) != std::string::npos;
	return has_prefix && has_mid;
}

bool is_a(std::string instance, value_ref type)
{
	return is_instance(instance)
		&& is_derived(
			variable(instance, "__class__")[0], list_cref(*list_ref(type)));
}

/*
Jam init code is:

# Copyright 2001, 2002, 2003 Dave Abrahams
# Copyright 2002, 2005 Rene Rivera
# Copyright 2002, 2003 Vladimir Prus
*/
const char * class_module::init_code = R"jam(

rule __test__ ( )
{
    import assert ;
    import "class" : new ;
    import errors : try catch ;

    # This will be the construction function for a class called 'myclass'.
    #
    class myclass
    {
        import assert ;

        rule __init__ ( x_ * : y_ * )
        {
            # Set some instance variables.
            x = $(x_) ;
            y = $(y_) ;
            foo += 10 ;
        }

        rule set-x ( newx * )
        {
            x = $(newx) ;
        }

        rule get-x ( )
        {
            return $(x) ;
        }

        rule set-y ( newy * )
        {
            y = $(newy) ;
        }

        rule get-y ( )
        {
            return $(y) ;
        }

        rule f ( )
        {
            return [ g $(x) ] ;
        }

        rule g ( args * )
        {
            if $(x) in $(y)
            {
                return $(x) ;
            }
            else if $(y) in $(x)
            {
                return $(y) ;
            }
            else
            {
                return ;
            }
        }

        rule get-class ( )
        {
            return $(__class__) ;
        }

        rule get-instance ( )
        {
            return $(__name__) ;
        }

        rule invariant ( )
        {
            assert.equal 1 : 1 ;
        }

        rule get-foo ( )
        {
            return $(foo) ;
        }
    }  # class myclass ;

    class derived1 : myclass
    {
        rule __init__ ( z_ )
        {
            myclass.__init__ $(z_) : X ;
            z = $(z_) ;
        }

        # Override g.
        #
        rule g ( args * )
        {
            return derived1.g ;
        }

        rule h ( )
        {
            return derived1.h ;
        }

        rule get-z ( )
        {
            return $(z) ;
        }

        # Check that 'assert.equal' visible in base class is visible here.
        #
        rule invariant2 ( )
        {
            assert.equal 2 : 2 ;
        }

        # Check that 'assert.variable-not-empty' visible in base class is
        # visible here.
        #
        rule invariant3 ( )
        {
            local v = 10 ;
            assert.variable-not-empty v ;
        }
    }  # class derived1 : myclass ;

    class derived2 : myclass
    {
        rule __init__ ( )
        {
            myclass.__init__ 1 : 2 ;
        }

        # Override g.
        #
        rule g ( args * )
        {
            return derived2.g ;
        }

        # Test the ability to call base class functions with qualification.
        #
        rule get-x ( )
        {
            return [ myclass.get-x ] ;
        }
    }  # class derived2 : myclass ;

    class derived2a : derived2
    {
        rule __init__
        {
            derived2.__init__ ;
        }
    }  # class derived2a : derived2 ;

    local rule expect_derived2 ( [derived2] x ) { }

    local a = [ new myclass 3 4 5 : 4 5 ] ;
    local b = [ new derived1 4 ] ;
    local b2 = [ new derived1 4 ] ;
    local c = [ new derived2 ] ;
    local d = [ new derived2 ] ;
    local e = [ new derived2a ] ;

    expect_derived2 $(d) ;
    expect_derived2 $(e) ;

    # Argument checking is set up to call exit(1) directly on failure, and we
    # can not hijack that with try, so we should better not do this test by
    # default. We could fix this by having errors look up and invoke the EXIT
    # rule instead; EXIT can be hijacked (;-)
    if --fail-typecheck in [ modules.peek : ARGV ]
    {
        try ;
        {
            expect_derived2 $(a) ;
        }
        catch
            "Expected an instance of derived2 but got" instead
            ;
    }

    #try ;
    #{
    #    new bad_subclass ;
    #}
    #catch
    #    bad_subclass.bad_subclass failed to call base class constructor
    #        myclass.__init__
    #  ;

    #try ;
    #{
    #    class bad_subclass ;
    #}
    #catch bad_subclass has already been declared ;

    assert.result 3 4 5 : $(a).get-x ;
    assert.result 4 5 : $(a).get-y ;
    assert.result 4 : $(b).get-x ;
    assert.result X : $(b).get-y ;
    assert.result 4 : $(b).get-z ;
    assert.result 1 : $(c).get-x ;
    assert.result 2 : $(c).get-y ;
    assert.result 4 5 : $(a).f ;
    assert.result derived1.g : $(b).f ;
    assert.result derived2.g : $(c).f ;
    assert.result derived2.g : $(d).f ;

    assert.result 10 : $(b).get-foo ;

    $(a).invariant ;
    $(b).invariant2 ;
    $(b).invariant3 ;

    # Check that the __class__ attribute is getting properly set.
    assert.result myclass : $(a).get-class ;
    assert.result derived1 : $(b).get-class ;
    assert.result $(a) : $(a).get-instance ;

    $(a).set-x a.x ;
    $(b).set-x b.x ;
    $(c).set-x c.x ;
    $(d).set-x d.x ;
    assert.result a.x : $(a).get-x ;
    assert.result b.x : $(b).get-x ;
    assert.result c.x : $(c).get-x ;
    assert.result d.x : $(d).get-x ;

    class derived3 : derived1 derived2
    {
        rule __init__ ( )
        {
        }
    }

    assert.result : bases myclass ;
    assert.result myclass : bases derived1 ;
    assert.result myclass : bases derived2 ;
    assert.result derived1 derived2 : bases derived3 ;

    assert.true is-derived derived1 : myclass ;
    assert.true is-derived derived2 : myclass ;
    assert.true is-derived derived3 : derived1 ;
    assert.true is-derived derived3 : derived2 ;
    assert.true is-derived derived3 : derived1 derived2 myclass ;
    assert.true is-derived derived3 : myclass ;

    assert.false is-derived myclass : derived1 ;

    assert.true is-instance $(a) ;
    assert.false is-instance bar ;

    assert.true is-a $(a) : myclass ;
    assert.true is-a $(c) : derived2 ;
    assert.true is-a $(d) : myclass ;
    assert.false is-a literal : myclass ;
}

)jam";

}}} // namespace b2::jam::klass
