#Called from the Galaxy Tool XML file
#import sys

def validate_input(trans, error_map, param_values, page_param_map):
    """Validates the min_size/max_size user input, before execution."""
    err_list = []
    for i, read_group in enumerate(param_values["read_group"]):
        err = dict()
        segments = read_group["segments"]
        if str(segments["type"]) != "paired":
            err_list.append(dict())
            continue

        min_size = str(segments["min_size"]).strip()
        max_size = str(segments["max_size"]).strip()
        #sys.stderr.write("DEBUG min_size=%r, max_size=%r\n" % (min_size, max_size))

        if min_size=="" and max_size=="":
            #Both missing is good
            pass
        elif min_size=="":
            err["min_size"] = "Minimum size required if maximum size given"
        elif max_size=="":
            err["max_size"] = "Maximum size required if minimum size given"
            
        if min_size:
            try:
                min_size_int = int(min_size)
                if min_size_int < 0:
                    err["min_size"] = "Minumum size must not be negative (%i)" % min_size_int
                    min_size = None # Avoid doing comparison below
            except ValueError:
                err["min_size"] = "Minimum size is not an integer (%s)" % min_size
                min_size = None # Avoid doing comparison below

        if max_size:
            try:
                max_size_int = int(max_size)
                if max_size_int< 0:
                    err["max_size"] = "Maximum size must not be negative (%i)" % max_size_int
                    max_size = None # Avoid doing comparison below
            except ValueError:
                err["max_size"] = "Maximum size is not an integer (%s)" % max_size
                max_size = None # Avoid doing comparison below

        if min_size and max_size and min_size_int > max_size_int:
            msg = "Minimum size must be less than maximum size (%i vs %i)" % (min_size_int, max_size_int)
            err["min_size"] = msg
            err["max_size"] = msg

        if err:
            err_list.append({"segments":err})
        else:
            err_list.append(dict())

    if any(err_list):
        #Return an error map only if any readgroup gave errors
        error_map["read_group"] = err_list
