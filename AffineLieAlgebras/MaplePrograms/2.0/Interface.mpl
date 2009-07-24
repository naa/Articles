########################################################################
########################################################################
###
###	Filename: Interface.mpl
###
###	Copyright (c) Maplesoft, a division of Waterloo Maple Inc. 1999-2004
###
###	You are permitted to copy, modify and distribute this code, as long as this
###	copyright notice is prominently included and remains intact. If any
###	modifications were done, a prominent notice of the fact that the code has been
###	modified, as well as a list of the modifications, must also be included. To the
###	maximum extent permitted by applicable laws, this material is provided "as is"
###	without any warranty or condition of any kind.
###
########################################################################
########################################################################
###
###	Description: This file defines a small package `Interface' for
###	manipulating interface types by name. The package exports the
###	following procedures
###
###		define -- define a new interface
###		extend -- extend an existing interface
###		extends -- test whether one interface extends another
###		equivalent -- test whether two interfaces are equivalent
###		savelib -- save an interface to a repository
###		satisfies -- test whether a module satisfies an interface
###
###	This package also serves as an example of the use of the `load='
###	option for modules. In this example, the option is used to
###	trigger the calling of the local thunk `setup' at load time.
###	This procedures simply ensures that the global type `interface'
###	is installed.
###
########################################################################
########################################################################

#$include	<assert.mi>

Interface := module()
	description "a package for manipulating interfaces";
	global	`type/interface`;
	export	define,		# define an interface
		extend,		# extend an interface
		extends,	# test for an extension
		equivalent,	# test equivalence
		savelib,	# save an interface
		satisfies;	# test whether a module satisfies an interface
	local	gassign,	# assign to a global variable
		totype,		# convert from interface name to type
		toset,		# convert from interface name to a set
		setup;		# install `type/interface` globally
	option	load = setup, package;

	# Define a global type for interfaces.
	# This assignment takes care of installing the type
	# in the Maple session in which this module definition
	# is evaluated. Calling `setup()' ensures that this type
	# installation also happens when the instantiated module
	# is read from a repository.
	`type/interface` := 'specfunc( { symbol, `::` }, `module` )';

	# Ensure that `type/interface` is defined. This thunk is
	# called when the instantiated `Interface' module is read
	# from a Maple repository.
	setup := proc()
		global `type/interface`;
		`type/interface` := 'specfunc( { symbol, `::` }, `module` )';
		NULL # quiet return
	end proc;

	# Assign to the global instance of a name
	gassign := proc( nom::symbol, val )
		option inline;
		eval( subs( _X = nom, proc() global _X; _X := val end ) )()
	end proc;

	# Convert an interface name to the corresponding type.
	totype := ( ifc::symbol ) -> ( `type/` || ifc );

	# Convert an interface name to a set of symbols (its exports).
	toset := ( ifc::symbol ) -> { op( ( `type/` || ifc ) ) };

	# Install a new interface into the type system.
	define := proc( ifc )
		description "define an interface";
		if not andmap( type, { args }, 'symbol' ) then
			error "arguments must all be symbols"
		end if;
		gassign( `type/` || ifc,
			 '`module`'( args[ 2 .. nargs ] ) );
		ifc # return the interface name
	end proc;

	# Implement subtyping.
	extend := proc( new, old )
		description "extend an existing interface";
		if not andmap( type, { args }, 'symbol' ) then
			error "arguments must all be symbols"
		end if;
		if not type( totype( old ), 'interface' ) then
			error "cannot find an interface named %1", old
		end if;
		define( new, op( totype( old ) ), args[3..nargs] )
	end proc;

	# Test whether ifc2 is an extension of ifc1.
	extends := proc( ifc1, ifc2 )
		description "test whether the second interface extends the first";
		local t1, t2;
		t1, t2 := op( map( totype, [ ifc1, ifc2 ] ) );
		if not type( [t1,t2], '[interface,interface]' ) then
			if not type( t1, 'interface' ) then
				error "arguments must be interface names, but got %1", ifc1
			else
				error "arguments must be interface names, but got %1", ifc1
			end if
		end if;
		toset( ifc1 ) subset toset( ifc2 )
	end proc;

	# Save an interface to the repository.
	savelib := proc()
		description "save a named interface to a repository";
		local ifc;
		for ifc in map( totype, [ args ] ) do
			if not type( ifc, 'interface' ) then
				error "arguments must be interfaces, but got %1", ifc
			end if;
			:-savelib( totype( ifc ) )
		end do
	end proc;

	# Test whether a module satisfies an interface. This is simply
	# an alternative to a call to `type()'.
	satisfies := proc( m, ifc )
		description "test whether a module satisfies an interface";
		if not type( totype( ifc ), 'interface' ) then
			error "second argument must be an interface name, "
				"but got %1", ifc
		end if;
		type( m, ifc )
	end proc;

	# Test whether two interfaces are equivalent. Since unevaluated
	# function calls compare differently if their arguments are
	# in a different order, we convert them to sets first, and
	# then test for equality.
	equivalent := proc( ifc1, ifc2 )
		description "test whether two interfaces are equivalent";
		local	t1, t2;
		t1, t2 := totype( ifc1 ), totype( ifc2 );
		if not type( t1, 'interface' ) then
			error "expecting an interface name, but got %1", ifc1
		elif not type( t2, 'interface' ) then
			error "expecting an interface name, but got %1", ifc2
		end if;
		evalb( { op( t1 ) } = { op( t2 ) } )
	end proc;
end module:
savelib( 'Interface' ):
#$ifdef	LOAD_REPOSITORY
#savelib( 'Interface' ):
#$endif	# LOAD_REPOSITORY

########################################################################
########################################################################
###
###	End of File: Interface.mpl
###
########################################################################
########################################################################
