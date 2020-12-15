import inspect
import os, re, stat
import sys
import simenv

import time

# dprint constants

global log
log = None

USE_VERBOSE = 0
ALWAYS_PRINT = 1
LOG_ONLY = 2

# There are three levels of error and warning output: errors are shown
# all the time and abort the code; warnings are shown all the time;
# informational messages ("debug messages") are only shown when
# --verbose is given.

def dprint(statement, flag=0, stderr=False, noNewLine=False):
    global log
     
    log_date = time.strftime("%Y-%m-%d %H:%M:%S")

    if not simenv.INITIALIZED:
        if stderr:
          dest = sys.stderr
        else:
          dest = sys.stdout
        dest.write("%s" % statement)
        if not noNewLine:
            dest.write("\n")
        return
        
    if log == None:
        log_dir = simenv.LOG_PATH
        
        if not(os.path.exists(log_dir)):
            try:
                os.makedirs(log_dir)
            except OSError, e:
                print "could not create logging directory \"%s\", %s" % (log_dir, e)
                sys.exit(1)
        
        log_path = os.sep.join([log_dir, "simfactory.log"])
        
        try:
            log = open(log_path, "a+")
        except:
            print "Could not open log file \"%s\" for writing" % log_path
            # Do not abort when the log file is not readable -- output
            # to stderr instead
            #sys.exit(1)
            log = sys.stderr
    
    
    log.write("[LOG:%s] %s\n" % (log_date, statement))
    # TODO: don't flush for every message; this is very slow
    log.flush()
    
    if flag == LOG_ONLY:
        return
    
    if flag == ALWAYS_PRINT or simenv.VERBOSE:
    
        if stderr:
            sys.stdout.flush()
            #stderr.write() doesn't append newline like print does.
            if noNewLine:
                sys.stderr.write("%s" % statement)
            else:
                sys.stderr.write("%s\n" % statement)
        else:
            if noNewLine:
                print statement,
            else:
                print statement
                
def fatal(message):
    global ALWAYS_PRINT
    global log
    
    message = "Error: %s" % message
    dprint(message, ALWAYS_PRINT, stderr=True)
    dprint("Aborting Simfactory.", ALWAYS_PRINT, stderr=True)
    
    if log != None:
        log.close()
        
    sys.exit(1)
    
def warning(message, noNewLine=False):
    global ALWAYS_PRINT
    
    message = "Warning: %s" % message
    dprint(message, ALWAYS_PRINT, stderr=True, noNewLine=noNewLine)
    sys.stderr.flush()

def info(message, noNewLine=False):
    
    message = "Info: %s" % message
    dprint(message, USE_VERBOSE, stderr=False, noNewLine=noNewLine)

def logonly(message, noNewLine=False):
    global LOG_ONLY

    dprint(message, LOG_ONLY, stderr=False, noNewLine=noNewLine)

def display(message, noNewLine=False):
    global ALWAYS_PRINT
    
    dprint(message, ALWAYS_PRINT, stderr=False, noNewLine=noNewLine)

        
def LineNumber():
    return inspect.currentframe().f_back.f_lineno
    
def FileExists(file):
    if file == None:
        return False
    return os.path.exists(file)
    
def MakeFileExecutable(file):
    perms = os.stat(file).st_mode
    perms |= stat.S_IRWXU
    if perms | stat.S_IRGRP:
        perms |= stat.S_IXGRP
    if perms | stat.S_IROTH:
        perms |= stat.S_IXOTH
    os.chmod(file, perms)
    

def Matches(pattern, value):
    p = re.compile(pattern)
    m = p.search(value)
    
    return m != None

def BuildWithSpace(str, numspaces):
    num_spaces = numspaces - len(str)
    
    fullString = str
    
    while num_spaces > 0:
        fullString = "%s " % fullString
        num_spaces = num_spaces - 1
        
    return fullString

def BuildStartEndString(start, end, numspaces):
    num_spaces = numspaces - len(start)
    
    fullString = start
    
    while num_spaces > 0:
        fullString = "%s " % fullString
        num_spaces = num_spaces - 1
        
    return "%s%s" % (fullString, end)
    
def PrintManyLeadingSpace(str, leading):
    lines = str.split("\n")
    
    for line in lines:
        sp = ""
        while leading > 0:
            sp = " %s" % sp
            leading = leading -1

        dprint("%s%s" % (sp, line))

# This function is used only in very few places; it should probably be
# removed, or moved to simopts. Furthermore, it should probably only
# translate strings to boolean values, and it is unclear whether the
# empty string should represent false. Maybe it should also know about
# true_values, and report an error when an unknown value is
# encountered.
def CoerceBool(val):
        false_values = ['False', 'false', 'No', 'no', 'F', 'f', '0', 0, False,
                        None, '']
        return val not in false_values
