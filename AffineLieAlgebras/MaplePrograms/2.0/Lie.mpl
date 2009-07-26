libname:="./lib",libname;

with(Interface);

define('Vectorspace'
       ,
       'dimension'
      );
AbstractVectorspace:=proc(dim)
    module()
    export dimension;
        dimension:=dim;
    end;
end;

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
       'representation',
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
define('Representation'
       ,
       'weights',
       'algebra'
       );

extend('HighestWeightRepresentation','Representation'
       ,
       'highestWeight',
       'mainChamber'
      );




anomalous_branching_coefficient:=proc(lambda::weight,
                                      Lmu::HighestWeightModule)
                    local i,j,k;
end;


