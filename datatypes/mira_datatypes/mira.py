"""
MiraAssemblyFormat class for the 'mira' format within Galaxy
"""

from galaxy.datatypes.data import Text


class MiraAssemblyFormat(Text):
    """MIRA Assembly Format  data"""
    file_ext = "mira"

    def sniff( self, filename ):
        """Determines whether the file is a MIRA Assembly Format file.

        Note currently this only detects MIRA Assembly Format v2.0,
        as used in MIRA v3.9 and v4.0.

        It does not detect MIRA Assembly Format v1 as used in both
        MIRA v3.2 and v3.4.
        """
        h = open(filename)
        line = h.readline()
        if line.rstrip() != "@Version\t2\t0":
            h.close()
            return False
        line = h.readline()
        if line.rstrip() != "@Program\tMIRALIB":
            h.close()
            return False
        return True
    
    def merge(split_files, output_file):
        """Merging multiple MIRA files is non-trivial and may not be possible..."""
        if len(split_files) == 1:
            #For one file only, use base class method (move/copy)
            return Text.merge(split_files, output_file)
        if not split_files:
            raise ValueError("Given no MIRA, %r, to merge into %s" \
                             % (split_files, output_file))
        raise NotImplementedError("Merging MIRA Assembly Files has not been implemented")
    merge = staticmethod(merge)

    def split( cls, input_datasets, subdir_generator_function, split_params):
        """Split a MIRA Assembly File (not implemented for now)."""
        if split_params is None:
            return None
        raise NotImplementedError("Can't yet split a MIRA Assembly Format file")
    merge = staticmethod(merge)
