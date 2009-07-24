libname:="./lib",libname;

with(Interface);

define('Vectorspace'
       ,
       'dimension'
      );

define('Basis'
       ,
       'coeffs',
       'vectorspace',
       'unitnames'
      );
define('Vector'
       ,
       'coordinates',
       'basis'
      );
extend('Weight','Vector'
       ,
       'module',
       'multilicity'
       );
extend('Root','Weight');
define('WeylReflection'
       ,
       'act'
       );
define('WeylGroup'
       ,
       'elements'
       );
define('Algebra'
       ,
       'weylGroup',
       'simpleRoots'
      );
extend('SubAlgebra','Algebra'
       ,
       'superGroup'
      );
define('Module'
       ,
       'weights',
       'algebra'
       );

extend('HighestWeightModule','Module'
       ,
       'highestWeight',
       'mainChamber'
      );

anomalous_branching_coefficient:=proc(lambda::weight,
                                      Lmu::HighestWeightModule)
                    local i,j,k;
end;


