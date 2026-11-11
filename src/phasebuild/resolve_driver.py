"""Isolated RESOLVE API driver for an explicit authorized bundled or external runtime."""
from __future__ import print_function
import sys
from solve_resolve import resolve
from solve_resolve.resolve import helper
from boost_adaptbx.boost.python import streambuf
from six.moves import cStringIO as StringIO

if __name__=='__main__':
    errors=StringIO()
    result=resolve.call(sys.stdin.read(),err_streambuf=streambuf(errors,buffer_size=0),
        out_streambuf=streambuf(sys.stdout,buffer_size=0),isizeit=12,n_crystal_map_max=1,resolve_libs=helper.resolve_libs())
    print('RESOLVE_API_STATUS',result.exit_info.status)
    if errors.getvalue():print(errors.getvalue())
    sys.exit(0 if result.exit_info.status==0 else 1)
