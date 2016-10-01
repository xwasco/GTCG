BENSOLVE is a MATLAB implementation of Benson's algorithm to
solve linear vector optimization problems of the form

   C-minimize Px
   s.t. Bx>=b

INSTALLATION: 

1. Unpack 'bensolve-x.x.zip'.
2. Run MatLab and change to directory 'bensolve-x.x'.
3. Run 'help bensolve' for further instructions.

SYSTEM REQUIREMENTS:

This program has been tested with MatLab 7.13.0.564 (R2011b). For several
(but not all) systems it works with precompiled mex-files of LP solvers
'cdd' and 'glpk', which are included in the corresponding subfolders.

LICENSE:

Copyright (C) 2011-2012  Andreas Löhne

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see http://www.gnu.org/licenses/.