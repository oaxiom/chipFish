"""

statistics.py

generic statistics implementations and other pass-through goodies from scipy.

Part of glbase.

"""

import config, utils

# Copyright (c) 1999-2007 Gary Strangman; All Rights Reserved.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# Comments and/or additions are welcome (send e-mail to:
# strang@nmr.mgh.harvard.edu).
#
"""
pstat.py module

#################################################
#######  Written by:  Gary Strangman  ###########
#######  Last modified:  Dec 18, 2007 ###########
#################################################

This module provides some useful list and array manipulation routines
modeled after those found in the |Stat package by Gary Perlman, plus a
number of other useful list/file manipulation functions.  The list-based
functions include:

      abut (source,*args)
      simpleabut (source, addon)
      colex (listoflists,cnums)
      collapse (listoflists,keepcols,collapsecols,fcn1=None,fcn2=None,cfcn=None)
      dm (listoflists,criterion)
      flat (l)
      linexand (listoflists,columnlist,valuelist)
      linexor (listoflists,columnlist,valuelist)
      linedelimited (inlist,delimiter)
      lineincols (inlist,colsize)
      lineincustcols (inlist,colsizes)
      list2string (inlist)
      makelol(inlist)
      makestr(x)
      printcc (lst,extra=2)
      printincols (listoflists,colsize)
      pl (listoflists)
      printl(listoflists)
      replace (lst,oldval,newval)
      recode (inlist,listmap,cols='all')
      remap (listoflists,criterion)
      roundlist (inlist,num_digits_to_round_floats_to)
      sortby(listoflists,sortcols)
      unique (inlist)
      duplicates(inlist)
      writedelimited (listoflists, delimiter, file, writetype='w')

Some of these functions have alternate versions which are defined only if
Numeric (NumPy) can be imported.  These functions are generally named as
above, with an 'a' prefix.

      aabut (source, *args)
      acolex (a,indices,axis=1)
      acollapse (a,keepcols,collapsecols,sterr=0,ns=0)
      adm (a,criterion)
      alinexand (a,columnlist,valuelist)
      alinexor (a,columnlist,valuelist)
      areplace (a,oldval,newval)
      arecode (a,listmap,col='all')
      arowcompare (row1, row2)
      arowsame (row1, row2)
      asortrows(a,axis=0)
      aunique(inarray)
      aduplicates(inarray)

Currently, the code is all but completely un-optimized.  In many cases, the
array versions of functions amount simply to aliases to built-in array
functions/methods.  Their inclusion here is for function name consistency.
"""

## CHANGE LOG:
## ==========
## 07-11-26 ... edited to work with numpy
## 01-11-15 ... changed list2string() to accept a delimiter
## 01-06-29 ... converted exec()'s to eval()'s to make compatible with Py2.1
## 01-05-31 ... added duplicates() and aduplicates() functions
## 00-12-28 ... license made GPL, docstring and import requirements
## 99-11-01 ... changed version to 0.3
## 99-08-30 ... removed get, getstrings, put, aget, aput (into io.py)
## 03/27/99 ... added areplace function, made replace fcn recursive
## 12/31/98 ... added writefc function for ouput to fixed column sizes
## 12/07/98 ... fixed import problem (failed on collapse() fcn)
##              added __version__ variable (now 0.2)
## 12/05/98 ... updated doc-strings
##              added features to collapse() function
##              added flat() function for lists
##              fixed a broken asortrows()
## 11/16/98 ... fixed minor bug in aput for 1D arrays
##
## 11/08/98 ... fixed aput to output large arrays correctly

import stats  # required 3rd party module
import string, copy
from types import *

__version__ = 0.4

###===========================  LIST FUNCTIONS  ==========================
###
### Here are the list functions, DEFINED FOR ALL SYSTEMS.
### Array functions (for NumPy-enabled computers) appear below.
###

def abut (source,*args):
    """
Like the |Stat abut command.  It concatenates two lists side-by-side
and returns the result.  '2D' lists are also accomodated for either argument
(source or addon).  CAUTION:  If one list is shorter, it will be repeated
until it is as long as the longest list.  If this behavior is not desired,
use pstat.simpleabut().

Usage:   abut(source, args)   where args=any # of lists
Returns: a list of lists as long as the LONGEST list past, source on the
         'left', lists in <args> attached consecutively on the 'right'
"""

    if type(source) not in [ListType,TupleType]:
        source = [source]
    for addon in args:
        if type(addon) not in [ListType,TupleType]:
            addon = [addon]
        if len(addon) < len(source):                # is source list longer?
            if len(source) % len(addon) == 0:        # are they integer multiples?
                repeats = len(source)/len(addon)    # repeat addon n times
                origadd = copy.deepcopy(addon)
                for i in range(repeats-1):
                    addon = addon + origadd
            else:
                repeats = len(source)/len(addon)+1  # repeat addon x times,
                origadd = copy.deepcopy(addon)      #    x is NOT an integer
                for i in range(repeats-1):
                    addon = addon + origadd
                    addon = addon[0:len(source)]
        elif len(source) < len(addon):                # is addon list longer?
            if len(addon) % len(source) == 0:        # are they integer multiples?
                repeats = len(addon)/len(source)    # repeat source n times
                origsour = copy.deepcopy(source)
                for i in range(repeats-1):
                    source = source + origsour
            else:
                repeats = len(addon)/len(source)+1  # repeat source x times,
                origsour = copy.deepcopy(source)    #   x is NOT an integer
                for i in range(repeats-1):
                    source = source + origsour
                source = source[0:len(addon)]

        source = simpleabut(source,addon)
    return source


def simpleabut (source, addon):
    """
Concatenates two lists as columns and returns the result.  '2D' lists
are also accomodated for either argument (source or addon).  This DOES NOT
repeat either list to make the 2 lists of equal length.  Beware of list pairs
with different lengths ... the resulting list will be the length of the
FIRST list passed.

Usage:   simpleabut(source,addon)  where source, addon=list (or list-of-lists)
Returns: a list of lists as long as source, with source on the 'left' and
                 addon on the 'right'
"""
    if type(source) not in [ListType,TupleType]:
        source = [source]
    if type(addon) not in [ListType,TupleType]:
        addon = [addon]
    minlen = min(len(source),len(addon))
    list = copy.deepcopy(source)                # start abut process
    if type(source[0]) not in [ListType,TupleType]:
        if type(addon[0]) not in [ListType,TupleType]:
            for i in range(minlen):
                list[i] = [source[i]] + [addon[i]]        # source/addon = column
        else:
            for i in range(minlen):
                list[i] = [source[i]] + addon[i]        # addon=list-of-lists
    else:
        if type(addon[0]) not in [ListType,TupleType]:
            for i in range(minlen):
                list[i] = source[i] + [addon[i]]        # source=list-of-lists
        else:
            for i in range(minlen):
                list[i] = source[i] + addon[i]        # source/addon = list-of-lists
    source = list
    return source


def colex (listoflists,cnums):
    """
Extracts from listoflists the columns specified in the list 'cnums'
(cnums can be an integer, a sequence of integers, or a string-expression that
corresponds to a slice operation on the variable x ... e.g., 'x[3:]' will colex
columns 3 onward from the listoflists).

Usage:   colex (listoflists,cnums)
Returns: a list-of-lists corresponding to the columns from listoflists
         specified by cnums, in the order the column numbers appear in cnums
"""
    global index
    column = 0
    if type(cnums) in [ListType,TupleType]:   # if multiple columns to get
        index = cnums[0]
        column = map(lambda x: x[index], listoflists)
        for col in cnums[1:]:
            index = col
            column = abut(column,map(lambda x: x[index], listoflists))
    elif type(cnums) == StringType:              # if an 'x[3:]' type expr.
        evalstring = 'map(lambda x: x'+cnums+', listoflists)'
        column = eval(evalstring)
    else:                                     # else it's just 1 col to get
        index = cnums
        column = map(lambda x: x[index], listoflists)
    return column


def collapse (listoflists,keepcols,collapsecols,fcn1=None,fcn2=None,cfcn=None):
     """
Averages data in collapsecol, keeping all unique items in keepcols
(using unique, which keeps unique LISTS of column numbers), retaining the
unique sets of values in keepcols, the mean for each.  Setting fcn1
and/or fcn2 to point to a function rather than None (e.g., stats.sterr, len)
will append those results (e.g., the sterr, N) after each calculated mean.
cfcn is the collapse function to apply (defaults to mean, defined here in the
pstat module to avoid circular imports with stats.py, but harmonicmean or
others could be passed).

Usage:    collapse (listoflists,keepcols,collapsecols,fcn1=None,fcn2=None,cfcn=None)
Returns: a list of lists with all unique permutations of entries appearing in
     columns ("conditions") specified by keepcols, abutted with the result of
     cfcn (if cfcn=None, defaults to the mean) of each column specified by
     collapsecols.
"""
     def collmean (inlist):
         s = 0
         for item in inlist:
             s = s + item
         return s/float(len(inlist))

     if type(keepcols) not in [ListType,TupleType]:
         keepcols = [keepcols]
     if type(collapsecols) not in [ListType,TupleType]:
         collapsecols = [collapsecols]
     if cfcn == None:
         cfcn = collmean
     if keepcols == []:
         means = [0]*len(collapsecols)
         for i in range(len(collapsecols)):
             avgcol = colex(listoflists,collapsecols[i])
             means[i] = cfcn(avgcol)
             if fcn1:
                 try:
                     test = fcn1(avgcol)
                 except:
                     test = 'N/A'
                     means[i] = [means[i], test]
             if fcn2:
                 try:
                     test = fcn2(avgcol)
                 except:
                     test = 'N/A'
                 try:
                     means[i] = means[i] + [len(avgcol)]
                 except TypeError:
                     means[i] = [means[i],len(avgcol)]
         return means
     else:
         values = colex(listoflists,keepcols)
         uniques = unique(values)
         uniques.sort()
         newlist = []
         if type(keepcols) not in [ListType,TupleType]:  keepcols = [keepcols]
         for item in uniques:
             if type(item) not in [ListType,TupleType]:  item =[item]
             tmprows = linexand(listoflists,keepcols,item)
             for col in collapsecols:
                 avgcol = colex(tmprows,col)
                 item.append(cfcn(avgcol))
                 if fcn1 <> None:
                     try:
                         test = fcn1(avgcol)
                     except:
                         test = 'N/A'
                     item.append(test)
                 if fcn2 <> None:
                     try:
                         test = fcn2(avgcol)
                     except:
                         test = 'N/A'
                     item.append(test)
                 newlist.append(item)
         return newlist


def dm (listoflists,criterion):
    """
Returns rows from the passed list of lists that meet the criteria in
the passed criterion expression (a string as a function of x; e.g., 'x[3]>=9'
will return all rows where the 4th column>=9 and "x[2]=='N'" will return rows
with column 2 equal to the string 'N').

Usage:   dm (listoflists, criterion)
Returns: rows from listoflists that meet the specified criterion.
"""
    function = 'filter(lambda x: '+criterion+',listoflists)'
    lines = eval(function)
    return lines


def flat(l):
    """
Returns the flattened version of a '2D' list.  List-correlate to the a.ravel()()
method of NumPy arrays.

Usage:    flat(l)
"""
    newl = []
    for i in range(len(l)):
        for j in range(len(l[i])):
            newl.append(l[i][j])
    return newl


def linexand (listoflists,columnlist,valuelist):
    """
Returns the rows of a list of lists where col (from columnlist) = val
(from valuelist) for EVERY pair of values (columnlist[i],valuelists[i]).
len(columnlist) must equal len(valuelist).

Usage:   linexand (listoflists,columnlist,valuelist)
Returns: the rows of listoflists where columnlist[i]=valuelist[i] for ALL i
"""
    if type(columnlist) not in [ListType,TupleType]:
        columnlist = [columnlist]
    if type(valuelist) not in [ListType,TupleType]:
        valuelist = [valuelist]
    criterion = ''
    for i in range(len(columnlist)):
        if type(valuelist[i])==StringType:
            critval = '\'' + valuelist[i] + '\''
        else:
            critval = str(valuelist[i])
        criterion = criterion + ' x['+str(columnlist[i])+']=='+critval+' and'
    criterion = criterion[0:-3]         # remove the "and" after the last crit
    function = 'filter(lambda x: '+criterion+',listoflists)'
    lines = eval(function)
    return lines


def linexor (listoflists,columnlist,valuelist):
    """
Returns the rows of a list of lists where col (from columnlist) = val
(from valuelist) for ANY pair of values (colunmlist[i],valuelist[i[).
One value is required for each column in columnlist.  If only one value
exists for columnlist but multiple values appear in valuelist, the
valuelist values are all assumed to pertain to the same column.

Usage:   linexor (listoflists,columnlist,valuelist)
Returns: the rows of listoflists where columnlist[i]=valuelist[i] for ANY i
"""
    if type(columnlist) not in [ListType,TupleType]:
        columnlist = [columnlist]
    if type(valuelist) not in [ListType,TupleType]:
        valuelist = [valuelist]
    criterion = ''
    if len(columnlist) == 1 and len(valuelist) > 1:
        columnlist = columnlist*len(valuelist)
    for i in range(len(columnlist)):          # build an exec string
        if type(valuelist[i])==StringType:
            critval = '\'' + valuelist[i] + '\''
        else:
            critval = str(valuelist[i])
        criterion = criterion + ' x['+str(columnlist[i])+']=='+critval+' or'
    criterion = criterion[0:-2]         # remove the "or" after the last crit
    function = 'filter(lambda x: '+criterion+',listoflists)'
    lines = eval(function)
    return lines


def linedelimited (inlist,delimiter):
    """
Returns a string composed of elements in inlist, with each element
separated by 'delimiter.'  Used by function writedelimited.  Use '\t'
for tab-delimiting.

Usage:   linedelimited (inlist,delimiter)
"""
    outstr = ''
    for item in inlist:
        if type(item) <> StringType:
            item = str(item)
        outstr = outstr + item + delimiter
    outstr = outstr[0:-1]
    return outstr


def lineincols (inlist,colsize):
    """
Returns a string composed of elements in inlist, with each element
right-aligned in columns of (fixed) colsize.

Usage:   lineincols (inlist,colsize)   where colsize is an integer
"""
    outstr = ''
    for item in inlist:
        if type(item) <> StringType:
            item = str(item)
        size = len(item)
        if size <= colsize:
            for i in range(colsize-size):
                outstr = outstr + ' '
            outstr = outstr + item
        else:
            outstr = outstr + item[0:colsize+1]
    return outstr


def lineincustcols (inlist,colsizes):
    """
Returns a string composed of elements in inlist, with each element
right-aligned in a column of width specified by a sequence colsizes.  The
length of colsizes must be greater than or equal to the number of columns
in inlist.

Usage:   lineincustcols (inlist,colsizes)
Returns: formatted string created from inlist
"""
    outstr = ''
    for i in range(len(inlist)):
        if type(inlist[i]) <> StringType:
            item = str(inlist[i])
        else:
            item = inlist[i]
        size = len(item)
        if size <= colsizes[i]:
            for j in range(colsizes[i]-size):
                outstr = outstr + ' '
            outstr = outstr + item
        else:
            outstr = outstr + item[0:colsizes[i]+1]
    return outstr


def list2string (inlist,delimit=' '):
    """
Converts a 1D list to a single long string for file output, using
the string.join function.

Usage:   list2string (inlist,delimit=' ')
Returns: the string created from inlist
"""
    stringlist = map(makestr,inlist)
    return string.join(stringlist,delimit)


def makelol(inlist):
    """
Converts a 1D list to a 2D list (i.e., a list-of-lists).  Useful when you
want to use put() to write a 1D list one item per line in the file.

Usage:   makelol(inlist)
Returns: if l = [1,2,'hi'] then returns [[1],[2],['hi']] etc.
"""
    x = []
    for item in inlist:
        x.append([item])
    return x


def makestr (x):
    if type(x) <> StringType:
        x = str(x)
    return x


def printcc (lst,extra=2):
    """
Prints a list of lists in columns, customized by the max size of items
within the columns (max size of items in col, plus 'extra' number of spaces).
Use 'dashes' or '\\n' in the list-of-lists to print dashes or blank lines,
respectively.

Usage:   printcc (lst,extra=2)
Returns: None
"""
    if type(lst[0]) not in [ListType,TupleType]:
        lst = [lst]
    rowstokill = []
    list2print = copy.deepcopy(lst)
    for i in range(len(lst)):
        if lst[i] == ['\n'] or lst[i]=='\n' or lst[i]=='dashes' or lst[i]=='' or lst[i]==['']:
            rowstokill = rowstokill + [i]
    rowstokill.reverse()   # delete blank rows from the end
    for row in rowstokill:
        del list2print[row]
    maxsize = [0]*len(list2print[0])
    for col in range(len(list2print[0])):
        items = colex(list2print,col)
        items = map(makestr,items)
        maxsize[col] = max(map(len,items)) + extra
    for row in lst:
        if row == ['\n'] or row == '\n' or row == '' or row == ['']:
            print
        elif row == ['dashes'] or row == 'dashes':
            dashes = [0]*len(maxsize)
            for j in range(len(maxsize)):
                dashes[j] = '-'*(maxsize[j]-2)
            print lineincustcols(dashes,maxsize)
        else:
            print lineincustcols(row,maxsize)
    return None


def printincols (listoflists,colsize):
    """
Prints a list of lists in columns of (fixed) colsize width, where
colsize is an integer.

Usage:   printincols (listoflists,colsize)
Returns: None
"""
    for row in listoflists:
        print lineincols(row,colsize)
    return None


def pl (listoflists):
    """
Prints a list of lists, 1 list (row) at a time.

Usage:   pl(listoflists)
Returns: None
"""
    for row in listoflists:
        if row[-1] == '\n':
            print row,
        else:
            print row
    return None


def printl(listoflists):
    """Alias for pl."""
    pl(listoflists)
    return


def replace (inlst,oldval,newval):
    """
Replaces all occurrences of 'oldval' with 'newval', recursively.

Usage:   replace (inlst,oldval,newval)
"""
    lst = inlst*1
    for i in range(len(lst)):
        if type(lst[i]) not in [ListType,TupleType]:
            if lst[i]==oldval: lst[i]=newval
        else:
            lst[i] = replace(lst[i],oldval,newval)
    return lst


def recode (inlist,listmap,cols=None):
    """
Changes the values in a list to a new set of values (useful when
you need to recode data from (e.g.) strings to numbers.  cols defaults
to None (meaning all columns are recoded).

Usage:   recode (inlist,listmap,cols=None)  cols=recode cols, listmap=2D list
Returns: inlist with the appropriate values replaced with new ones
"""
    lst = copy.deepcopy(inlist)
    if cols != None:
        if type(cols) not in [ListType,TupleType]:
            cols = [cols]
        for col in cols:
            for row in range(len(lst)):
                try:
                    idx = colex(listmap,0).index(lst[row][col])
                    lst[row][col] = listmap[idx][1]
                except ValueError:
                    pass
    else:
        for row in range(len(lst)):
            for col in range(len(lst)):
                try:
                    idx = colex(listmap,0).index(lst[row][col])
                    lst[row][col] = listmap[idx][1]
                except ValueError:
                    pass
    return lst


def remap (listoflists,criterion):
    """
Remaps values in a given column of a 2D list (listoflists).  This requires
a criterion as a function of 'x' so that the result of the following is
returned ... map(lambda x: 'criterion',listoflists).

Usage:   remap(listoflists,criterion)    criterion=string
Returns: remapped version of listoflists
"""
    function = 'map(lambda x: '+criterion+',listoflists)'
    lines = eval(function)
    return lines


def roundlist (inlist,digits):
    """
Goes through each element in a 1D or 2D inlist, and applies the following
function to all elements of FloatType ... round(element,digits).

Usage:   roundlist(inlist,digits)
Returns: list with rounded floats
"""
    if type(inlist[0]) in [IntType, FloatType]:
        inlist = [inlist]
    l = inlist*1
    for i in range(len(l)):
        for j in range(len(l[i])):
            if type(l[i][j])==FloatType:
                l[i][j] = round(l[i][j],digits)
    return l


def sortby(listoflists,sortcols):
    """
Sorts a list of lists on the column(s) specified in the sequence
sortcols.

Usage:   sortby(listoflists,sortcols)
Returns: sorted list, unchanged column ordering
"""
    newlist = abut(colex(listoflists,sortcols),listoflists)
    newlist.sort()
    try:
        numcols = len(sortcols)
    except TypeError:
        numcols = 1
    crit = '[' + str(numcols) + ':]'
    newlist = colex(newlist,crit)
    return newlist

def nonrepeats(inlist):
    """
Returns items that are NOT duplicated in the first dim of the passed list.

Usage:   nonrepeats (inlist)
"""
    nonrepeats = []
    for i in range(len(inlist)):
        if inlist.count(inlist[i]) == 1:
            nonrepeats.append(inlist[i])
    return nonrepeats


# Copyright (c) 1999-2007 Gary Strangman; All Rights Reserved.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# Comments and/or additions are welcome (send e-mail to:
# strang@nmr.mgh.harvard.edu).
#
"""
stats.py module

(Requires pstat.py module.)

#################################################
#######  Written by:  Gary Strangman  ###########
#######  Last modified:  Dec 18, 2007 ###########
#################################################

A collection of basic statistical functions for python.  The function
names appear below.

IMPORTANT:  There are really *3* sets of functions.  The first set has an 'l'
prefix, which can be used with list or tuple arguments.  The second set has
an 'a' prefix, which can accept NumPy array arguments.  These latter
functions are defined only when NumPy is available on the system.  The third
type has NO prefix (i.e., has the name that appears below).  Functions of
this set are members of a "Dispatch" class, c/o David Ascher.  This class
allows different functions to be called depending on the type of the passed
arguments.  Thus, stats.mean is a member of the Dispatch class and
stats.mean(range(20)) will call stats.lmean(range(20)) while
stats.mean(Numeric.arange(20)) will call stats.amean(Numeric.arange(20)).
This is a handy way to keep consistent function names when different
argument types require different functions to be called.  Having
implementated the Dispatch class, however, means that to get info on
a given function, you must use the REAL function name ... that is
"print stats.lmean.__doc__" or "print stats.amean.__doc__" work fine,
while "print stats.mean.__doc__" will print the doc for the Dispatch
class.  NUMPY FUNCTIONS ('a' prefix) generally have more argument options
but should otherwise be consistent with the corresponding list functions.

Disclaimers:  The function list is obviously incomplete and, worse, the
functions are not optimized.  All functions have been tested (some more
so than others), but they are far from bulletproof.  Thus, as with any
free software, no warranty or guarantee is expressed or implied. :-)  A
few extra functions that don't appear in the list below can be found by
interested treasure-hunters.  These functions don't necessarily have
both list and array versions but were deemed useful

CENTRAL TENDENCY:  geometricmean
                   harmonicmean
                   mean
                   median
                   medianscore
                   mode

MOMENTS:  moment
          variation
          skew
          kurtosis
          skewtest   (for Numpy arrays only)
          kurtosistest (for Numpy arrays only)
          normaltest (for Numpy arrays only)

ALTERED VERSIONS:  tmean  (for Numpy arrays only)
                   tvar   (for Numpy arrays only)
                   tmin   (for Numpy arrays only)
                   tmax   (for Numpy arrays only)
                   tstdev (for Numpy arrays only)
                   tsem   (for Numpy arrays only)
                   describe

FREQUENCY STATS:  itemfreq
                  scoreatpercentile
                  percentileofscore
                  histogram
                  cumfreq
                  relfreq

VARIABILITY:  obrientransform
              samplevar
              samplestdev
              signaltonoise (for Numpy arrays only)
              var
              stdev
              sterr
              sem
              z
              zs
              zmap (for Numpy arrays only)

TRIMMING FCNS:  threshold (for Numpy arrays only)
                trimboth
                trim1
                round (round all vals to 'n' decimals; Numpy only)

CORRELATION FCNS:  covariance  (for Numpy arrays only)
                   correlation (for Numpy arrays only)
                   paired
                   pearsonr
                   spearmanr
                   pointbiserialr
                   kendalltau
                   linregress

INFERENTIAL STATS:  ttest_1samp
                    ttest_ind
                    ttest_rel
                    chisquare
                    ks_2samp
                    mannwhitneyu
                    ranksums
                    wilcoxont
                    kruskalwallish
                    friedmanchisquare

PROBABILITY CALCS:  chisqprob
                    erfcc
                    zprob
                    ksprob
                    fprob
                    betacf
                    gammln
                    betai

ANOVA FUNCTIONS:  F_oneway
                  F_value

SUPPORT FUNCTIONS:  writecc
                    incr
                    sign  (for Numpy arrays only)
                    sum
                    cumsum
                    ss
                    summult
                    sumdiffsquared
                    square_of_sums
                    shellsort
                    rankdata
                    outputpairedstats
                    findwithin
"""

import pstat               # required 3rd party module
import math, string, copy  # required python modules
from types import *

__version__ = 0.6

############# DISPATCH CODE ##############


class Dispatch:
    """
    The Dispatch class, care of David Ascher, allows different functions to
    be called depending on the argument types.  This way, there can be one
    function name regardless of the argument type.  To access function doc
    in stats.py module, prefix the function with an 'l' or 'a' for list or
    array arguments, respectively.  That is, print stats.lmean.__doc__ or
    print stats.amean.__doc__ or whatever.
    """

    def __init__(self, *tuples):
        self._dispatch = {}
        for func, types in tuples:
            for t in types:
                if t in self._dispatch.keys():
                    raise ValueError, "can't have two dispatches on "+str(t)
                self._dispatch[t] = func
        self._types = self._dispatch.keys()

    def __call__(self, arg1, *args, **kw):
        if type(arg1) not in self._types:
            raise TypeError, "don't know how to dispatch %s arguments" %  type(arg1)
        return apply(self._dispatch[type(arg1)], (arg1,) + args, kw)

    def __doc__(self):
        """
        Pass through at least one of the __doc__ strings, and
        add a warning that the args may not be correct.
        """
        keys = self._dispatch
        return(self._dispatch[key[0]].__doc__)


##########################################################################
########################   LIST-BASED FUNCTIONS   ########################
##########################################################################

### Define these regardless

####################################
#######  CENTRAL TENDENCY  #########
####################################

def lgeometricmean (inlist):
    """
    Calculates the geometric mean of the values in the passed list.
    That is:  n-th root of (x1 * x2 * ... * xn).  Assumes a '1D' list.

    Usage:   lgeometricmean(inlist)
    """
    mult = 1.0
    one_over_n = 1.0/len(inlist)
    for item in inlist:
        mult = mult * pow(item,one_over_n)
    return mult


def lharmonicmean (inlist):
    """
Calculates the harmonic mean of the values in the passed list.
That is:  n / (1/x1 + 1/x2 + ... + 1/xn).  Assumes a '1D' list.

Usage:   lharmonicmean(inlist)
"""
    sum = 0
    for item in inlist:
        sum = sum + 1.0/item
    return len(inlist) / sum


def lmean (inlist):
    """
Returns the arithematic mean of the values in the passed list.
Assumes a '1D' list, but will function on the 1st dim of an array(!).

Usage:   lmean(inlist)
"""
    sum = 0
    for item in inlist:
        sum = sum + item
    return sum/float(len(inlist))


def lmedian (inlist,numbins=1000):
    """
Returns the computed median value of a list of numbers, given the
number of bins to use for the histogram (more bins brings the computed value
closer to the median score, default number of bins = 1000).  See G.W.
Heiman's Basic Stats (1st Edition), or CRC Probability & Statistics.

Usage:   lmedian (inlist, numbins=1000)
"""
    (hist, smallest, binsize, extras) = histogram(inlist,numbins,[min(inlist),max(inlist)]) # make histog
    cumhist = cumsum(hist)            # make cumulative histogram
    for i in range(len(cumhist)):       # get 1st(!) index holding 50%ile score
        if cumhist[i]>=len(inlist)/2.0:
            cfbin = i
            break
    LRL = smallest + binsize*cfbin      # get lower read limit of that bin
    cfbelow = cumhist[cfbin-1]
    freq = float(hist[cfbin])               # frequency IN the 50%ile bin
    median = LRL + ((len(inlist)/2.0 - cfbelow)/float(freq))*binsize  # median formula
    return median


def lmedianscore (inlist):
    """
Returns the 'middle' score of the passed list.  If there is an even
number of scores, the mean of the 2 middle scores is returned.

Usage:   lmedianscore(inlist)
"""

    newlist = copy.deepcopy(inlist)
    newlist.sort()
    if len(newlist) % 2 == 0:   # if even number of scores, average middle 2
        index = len(newlist)/2  # integer division correct
        median = float(newlist[index] + newlist[index-1]) /2
    else:
        index = len(newlist)/2  # int divsion gives mid value when count from 0
        median = newlist[index]
    return median


def lmode(inlist):
    """
Returns a list of the modal (most common) score(s) in the passed
list.  If there is more than one such score, all are returned.  The
bin-count for the mode(s) is also returned.

Usage:   lmode(inlist)
Returns: bin-count for mode(s), a list of modal value(s)
"""

    scores = pstat.unique(inlist)
    scores.sort()
    freq = []
    for item in scores:
        freq.append(inlist.count(item))
    maxfreq = max(freq)
    mode = []
    stillmore = 1
    while stillmore:
        try:
            indx = freq.index(maxfreq)
            mode.append(scores[indx])
            del freq[indx]
            del scores[indx]
        except ValueError:
            stillmore=0
    return maxfreq, mode


####################################
############  MOMENTS  #############
####################################

def lmoment(inlist,moment=1):
    """
Calculates the nth moment about the mean for a sample (defaults to
the 1st moment).  Used to calculate coefficients of skewness and kurtosis.

Usage:   lmoment(inlist,moment=1)
Returns: appropriate moment (r) from ... 1/n * SUM((inlist(i)-mean)**r)
"""
    if moment == 1:
        return 0.0
    else:
        mn = mean(inlist)
        n = len(inlist)
        s = 0
        for x in inlist:
            s = s + (x-mn)**moment
        return s/float(n)


def lvariation(inlist):
    """
Returns the coefficient of variation, as defined in CRC Standard
Probability and Statistics, p.6.

Usage:   lvariation(inlist)
"""
    return 100.0*samplestdev(inlist)/float(mean(inlist))


def lskew(inlist):
    """
Returns the skewness of a distribution, as defined in Numerical
Recipies (alternate defn in CRC Standard Probability and Statistics, p.6.)

Usage:   lskew(inlist)
"""
    return moment(inlist,3)/pow(moment(inlist,2),1.5)


def lkurtosis(inlist):
    """
Returns the kurtosis of a distribution, as defined in Numerical
Recipies (alternate defn in CRC Standard Probability and Statistics, p.6.)

Usage:   lkurtosis(inlist)
"""
    return moment(inlist,4)/pow(moment(inlist,2),2.0)


def ldescribe(inlist):
    """
Returns some descriptive statistics of the passed list (assumed to be 1D).

Usage:   ldescribe(inlist)
Returns: n, mean, standard deviation, skew, kurtosis
"""
    n = len(inlist)
    mm = (min(inlist),max(inlist))
    m = mean(inlist)
    sd = stdev(inlist)
    sk = skew(inlist)
    kurt = kurtosis(inlist)
    return n, mm, m, sd, sk, kurt


####################################
#######  FREQUENCY STATS  ##########
####################################

def litemfreq(inlist):
    """
Returns a list of pairs.  Each pair consists of one of the scores in inlist
and it's frequency count.  Assumes a 1D list is passed.

Usage:   litemfreq(inlist)
Returns: a 2D frequency table (col [0:n-1]=scores, col n=frequencies)
"""
    scores = pstat.unique(inlist)
    scores.sort()
    freq = []
    for item in scores:
        freq.append(inlist.count(item))
    return pstat.abut(scores, freq)


def lscoreatpercentile (inlist, percent):
    """
Returns the score at a given percentile relative to the distribution
given by inlist.

Usage:   lscoreatpercentile(inlist,percent)
"""
    if percent > 1:
        print "\nDividing percent>1 by 100 in lscoreatpercentile().\n"
        percent = percent / 100.0
    targetcf = percent*len(inlist)
    h, lrl, binsize, extras = histogram(inlist)
    cumhist = cumsum(copy.deepcopy(h))
    for i in range(len(cumhist)):
        if cumhist[i] >= targetcf:
            break
    score = binsize * ((targetcf - cumhist[i-1]) / float(h[i])) + (lrl+binsize*i)
    return score


def lpercentileofscore (inlist, score,histbins=10,defaultlimits=None):
    """
Returns the percentile value of a score relative to the distribution
given by inlist.  Formula depends on the values used to histogram the data(!).

Usage:   lpercentileofscore(inlist,score,histbins=10,defaultlimits=None)
"""

    h, lrl, binsize, extras = histogram(inlist,histbins,defaultlimits)
    cumhist = cumsum(copy.deepcopy(h))
    i = int((score - lrl)/float(binsize))
    pct = (cumhist[i-1]+((score-(lrl+binsize*i))/float(binsize))*h[i])/float(len(inlist)) * 100
    return pct


def lhistogram (inlist,numbins=10,defaultreallimits=None,printextras=0):
    """
Returns (i) a list of histogram bin counts, (ii) the smallest value
of the histogram binning, and (iii) the bin width (the last 2 are not
necessarily integers).  Default number of bins is 10.  If no sequence object
is given for defaultreallimits, the routine picks (usually non-pretty) bins
spanning all the numbers in the inlist.

Usage:   lhistogram (inlist, numbins=10, defaultreallimits=None,suppressoutput=0)
Returns: list of bin values, lowerreallimit, binsize, extrapoints
"""
    if (defaultreallimits <> None):
        if type(defaultreallimits) not in [ListType,TupleType] or len(defaultreallimits)==1: # only one limit given, assumed to be lower one & upper is calc'd
            lowerreallimit = defaultreallimits
            upperreallimit = 1.000001 * max(inlist)
        else: # assume both limits given
            lowerreallimit = defaultreallimits[0]
            upperreallimit = defaultreallimits[1]
        binsize = (upperreallimit-lowerreallimit)/float(numbins)
    else:    # no limits given for histogram, both must be calc'd
        estbinwidth=(max(inlist)-min(inlist))/float(numbins) +1e-6 #1=>cover all
        binsize = ((max(inlist)-min(inlist)+estbinwidth))/float(numbins)
        lowerreallimit = min(inlist) - binsize/2 #lower real limit,1st bin
    bins = [0]*(numbins)
    extrapoints = 0
    for num in inlist:
        try:
            if (num-lowerreallimit) < 0:
                extrapoints = extrapoints + 1
            else:
                bintoincrement = int((num-lowerreallimit)/float(binsize))
                bins[bintoincrement] = bins[bintoincrement] + 1
        except:
            extrapoints = extrapoints + 1
    if (extrapoints > 0 and printextras == 1):
        print '\nPoints outside given histogram range =',extrapoints
    return (bins, lowerreallimit, binsize, extrapoints)


def lcumfreq(inlist,numbins=10,defaultreallimits=None):
    """
Returns a cumulative frequency histogram, using the histogram function.

Usage:   lcumfreq(inlist,numbins=10,defaultreallimits=None)
Returns: list of cumfreq bin values, lowerreallimit, binsize, extrapoints
"""
    h,l,b,e = histogram(inlist,numbins,defaultreallimits)
    cumhist = cumsum(copy.deepcopy(h))
    return cumhist,l,b,e


def lrelfreq(inlist,numbins=10,defaultreallimits=None):
    """
Returns a relative frequency histogram, using the histogram function.

Usage:   lrelfreq(inlist,numbins=10,defaultreallimits=None)
Returns: list of cumfreq bin values, lowerreallimit, binsize, extrapoints
"""
    h,l,b,e = histogram(inlist,numbins,defaultreallimits)
    for i in range(len(h)):
        h[i] = h[i]/float(len(inlist))
    return h,l,b,e


####################################
#####  VARIABILITY FUNCTIONS  ######
####################################

def lobrientransform(*args):
    """
Computes a transform on input data (any number of columns).  Used to
test for homogeneity of variance prior to running one-way stats.  From
Maxwell and Delaney, p.112.

Usage:   lobrientransform(*args)
Returns: transformed data for use in an ANOVA
"""
    TINY = 1e-10
    k = len(args)
    n = [0.0]*k
    v = [0.0]*k
    m = [0.0]*k
    nargs = []
    for i in range(k):
        nargs.append(copy.deepcopy(args[i]))
        n[i] = float(len(nargs[i]))
        v[i] = var(nargs[i])
        m[i] = mean(nargs[i])
    for j in range(k):
        for i in range(n[j]):
            t1 = (n[j]-1.5)*n[j]*(nargs[j][i]-m[j])**2
            t2 = 0.5*v[j]*(n[j]-1.0)
            t3 = (n[j]-1.0)*(n[j]-2.0)
            nargs[j][i] = (t1-t2) / float(t3)
    check = 1
    for j in range(k):
        if v[j] - mean(nargs[j]) > TINY:
            check = 0
    if check <> 1:
        raise ValueError, 'Problem in obrientransform.'
    else:
        return nargs


def lsamplevar (inlist):
    """
Returns the variance of the values in the passed list using
N for the denominator (i.e., DESCRIBES the sample variance only).

Usage:   lsamplevar(inlist)
"""
    n = len(inlist)
    mn = mean(inlist)
    deviations = []
    for item in inlist:
        deviations.append(item-mn)
    return ss(deviations)/float(n)


def lsamplestdev (inlist):
    """
Returns the standard deviation of the values in the passed list using
N for the denominator (i.e., DESCRIBES the sample stdev only).

Usage:   lsamplestdev(inlist)
"""
    return math.sqrt(samplevar(inlist))


def lcov (x,y, keepdims=0):
    """
Returns the estimated covariance of the values in the passed
array (i.e., N-1).  Dimension can equal None (ravel array first), an
integer (the dimension over which to operate), or a sequence (operate
over multiple dimensions).  Set keepdims=1 to return an array with the
same number of dimensions as inarray.

Usage:   lcov(x,y,keepdims=0)
"""

    n = len(x)
    xmn = mean(x)
    ymn = mean(y)
    xdeviations = [0]*len(x)
    ydeviations = [0]*len(y)
    for i in range(len(x)):
        xdeviations[i] = x[i] - xmn
        ydeviations[i] = y[i] - ymn
    ss = 0.0
    for i in range(len(xdeviations)):
        ss = ss + xdeviations[i]*ydeviations[i]
    return ss/float(n-1)


def lvar (inlist):
    """
Returns the variance of the values in the passed list using N-1
for the denominator (i.e., for estimating population variance).

Usage:   lvar(inlist)
"""
    n = len(inlist)
    mn = mean(inlist)
    deviations = [0]*len(inlist)
    for i in range(len(inlist)):
        deviations[i] = inlist[i] - mn
    return ss(deviations)/float(n-1)


def lstdev (inlist):
    """
Returns the standard deviation of the values in the passed list
using N-1 in the denominator (i.e., to estimate population stdev).

Usage:   lstdev(inlist)
"""
    return math.sqrt(var(inlist))


def lsterr(inlist):
    """
Returns the standard error of the values in the passed list using N-1
in the denominator (i.e., to estimate population standard error).

Usage:   lsterr(inlist)
"""
    return stdev(inlist) / float(math.sqrt(len(inlist)))


def lsem (inlist):
    """
Returns the estimated standard error of the mean (sx-bar) of the
values in the passed list.  sem = stdev / sqrt(n)

Usage:   lsem(inlist)
"""
    sd = stdev(inlist)
    n = len(inlist)
    return sd/math.sqrt(n)


def lz (inlist, score):
    """
Returns the z-score for a given input score, given that score and the
list from which that score came.  Not appropriate for population calculations.

Usage:   lz(inlist, score)
"""
    z = (score-mean(inlist))/samplestdev(inlist)
    return z


def lzs (inlist):
    """
Returns a list of z-scores, one for each score in the passed list.

Usage:   lzs(inlist)
"""
    zscores = []
    for item in inlist:
        zscores.append(z(inlist,item))
    return zscores


####################################
#######  TRIMMING FUNCTIONS  #######
####################################

def ltrimboth (l,proportiontocut):
    """
Slices off the passed proportion of items from BOTH ends of the passed
list (i.e., with proportiontocut=0.1, slices 'leftmost' 10% AND 'rightmost'
10% of scores.  Assumes list is sorted by magnitude.  Slices off LESS if
proportion results in a non-integer slice index (i.e., conservatively
slices off proportiontocut).

Usage:   ltrimboth (l,proportiontocut)
Returns: trimmed version of list l
"""
    lowercut = int(proportiontocut*len(l))
    uppercut = len(l) - lowercut
    return l[lowercut:uppercut]


def ltrim1 (l,proportiontocut,tail='right'):
    """
Slices off the passed proportion of items from ONE end of the passed
list (i.e., if proportiontocut=0.1, slices off 'leftmost' or 'rightmost'
10% of scores).  Slices off LESS if proportion results in a non-integer
slice index (i.e., conservatively slices off proportiontocut).

Usage:   ltrim1 (l,proportiontocut,tail='right')  or set tail='left'
Returns: trimmed version of list l
"""
    if tail == 'right':
        lowercut = 0
        uppercut = len(l) - int(proportiontocut*len(l))
    elif tail == 'left':
        lowercut = int(proportiontocut*len(l))
        uppercut = len(l)
    return l[lowercut:uppercut]


####################################
#####  CORRELATION FUNCTIONS  ######
####################################

def lpaired(x,y):
    """
Interactively determines the type of data and then runs the
appropriated statistic for paired group data.

Usage:   lpaired(x,y)
Returns: appropriate statistic name, value, and probability
"""
    samples = ''
    while samples not in ['i','r','I','R','c','C']:
        print '\nIndependent or related samples, or correlation (i,r,c): ',
        samples = raw_input()

    if samples in ['i','I','r','R']:
        print '\nComparing variances ...',
# USE O'BRIEN'S TEST FOR HOMOGENEITY OF VARIANCE, Maxwell & delaney, p.112
        r = obrientransform(x,y)
        f,p = F_oneway(pstat.colex(r,0),pstat.colex(r,1))
        if p<0.05:
            vartype='unequal, p='+str(round(p,4))
        else:
            vartype='equal'
        print vartype
        if samples in ['i','I']:
            if vartype[0]=='e':
                t,p = ttest_ind(x,y,0)
                print '\nIndependent samples t-test:  ', round(t,4),round(p,4)
            else:
                if len(x)>20 or len(y)>20:
                    z,p = ranksums(x,y)
                    print '\nRank Sums test (NONparametric, n>20):  ', round(z,4),round(p,4)
                else:
                    u,p = mannwhitneyu(x,y)
                    print '\nMann-Whitney U-test (NONparametric, ns<20):  ', round(u,4),round(p,4)

        else:  # RELATED SAMPLES
            if vartype[0]=='e':
                t,p = ttest_rel(x,y,0)
                print '\nRelated samples t-test:  ', round(t,4),round(p,4)
            else:
                t,p = ranksums(x,y)
                print '\nWilcoxon T-test (NONparametric):  ', round(t,4),round(p,4)
    else:  # CORRELATION ANALYSIS
        corrtype = ''
        while corrtype not in ['c','C','r','R','d','D']:
            print '\nIs the data Continuous, Ranked, or Dichotomous (c,r,d): ',
            corrtype = raw_input()
        if corrtype in ['c','C']:
            m,b,r,p,see = linregress(x,y)
            print '\nLinear regression for continuous variables ...'
            lol = [['Slope','Intercept','r','Prob','SEestimate'],[round(m,4),round(b,4),round(r,4),round(p,4),round(see,4)]]
            pstat.printcc(lol)
        elif corrtype in ['r','R']:
            r,p = spearmanr(x,y)
            print '\nCorrelation for ranked variables ...'
            print "Spearman's r: ",round(r,4),round(p,4)
        else: # DICHOTOMOUS
            r,p = pointbiserialr(x,y)
            print '\nAssuming x contains a dichotomous variable ...'
            print 'Point Biserial r: ',round(r,4),round(p,4)
    print '\n\n'
    return None


def lpearsonr(x,y):
    """
Calculates a Pearson correlation coefficient and the associated
probability value.  Taken from Heiman's Basic Statistics for the Behav.
Sci (2nd), p.195.

Usage:   lpearsonr(x,y)   where x and y are equal-length lists
Returns: Pearson's r value, two-tailed p-value
"""
    TINY = 1.0e-30
    if len(x) <> len(y):
        raise ValueError, 'Input values not paired in pearsonr.  Aborting.'
    n = len(x)
    x = map(float,x)
    y = map(float,y)
    xmean = mean(x)
    ymean = mean(y)
    r_num = n*(summult(x,y)) - sum(x)*sum(y)
    r_den = math.sqrt((n*ss(x) - square_of_sums(x))*(n*ss(y)-square_of_sums(y)))
    r = (r_num / r_den)  # denominator already a float
    df = n-2
    t = r*math.sqrt(df/((1.0-r+TINY)*(1.0+r+TINY)))
    prob = betai(0.5*df,0.5,df/float(df+t*t))
    return r, prob


def llincc(x,y):
    """
Calculates Lin's concordance correlation coefficient.

Usage:   alincc(x,y)    where x, y are equal-length arrays
Returns: Lin's CC
"""
    covar = lcov(x,y)*(len(x)-1)/float(len(x))  # correct denom to n
    xvar = lvar(x)*(len(x)-1)/float(len(x))  # correct denom to n
    yvar = lvar(y)*(len(y)-1)/float(len(y))  # correct denom to n
    lincc = (2 * covar) / ((xvar+yvar) +((amean(x)-amean(y))**2))
    return lincc


def lspearmanr(x,y):
    """
Calculates a Spearman rank-order correlation coefficient.  Taken
from Heiman's Basic Statistics for the Behav. Sci (1st), p.192.

Usage:   lspearmanr(x,y)      where x and y are equal-length lists
Returns: Spearman's r, two-tailed p-value
"""
    TINY = 1e-30
    if len(x) <> len(y):
        raise ValueError, 'Input values not paired in spearmanr.  Aborting.'
    n = len(x)
    rankx = rankdata(x)
    ranky = rankdata(y)
    dsq = sumdiffsquared(rankx,ranky)
    rs = 1 - 6*dsq / float(n*(n**2-1))
    t = rs * math.sqrt((n-2) / ((rs+1.0)*(1.0-rs)))
    df = n-2
    probrs = betai(0.5*df,0.5,df/(df+t*t))  # t already a float
# probability values for rs are from part 2 of the spearman function in
# Numerical Recipies, p.510.  They are close to tables, but not exact. (?)
    return rs, probrs


def lpointbiserialr(x,y):
    """
Calculates a point-biserial correlation coefficient and the associated
probability value.  Taken from Heiman's Basic Statistics for the Behav.
Sci (1st), p.194.

Usage:   lpointbiserialr(x,y)     where x,y are equal-length lists
Returns: Point-biserial r, two-tailed p-value
"""
    TINY = 1e-30
    if len(x) <> len(y):
        raise ValueError, 'INPUT VALUES NOT PAIRED IN pointbiserialr.  ABORTING.'
    data = pstat.abut(x,y)
    categories = pstat.unique(x)
    if len(categories) <> 2:
        raise ValueError, "Exactly 2 categories required for pointbiserialr()."
    else:   # there are 2 categories, continue
        codemap = pstat.abut(categories,range(2))
        recoded = pstat.recode(data,codemap,0)
        x = pstat.linexand(data,0,categories[0])
        y = pstat.linexand(data,0,categories[1])
        xmean = mean(pstat.colex(x,1))
        ymean = mean(pstat.colex(y,1))
        n = len(data)
        adjust = math.sqrt((len(x)/float(n))*(len(y)/float(n)))
        rpb = (ymean - xmean)/samplestdev(pstat.colex(data,1))*adjust
        df = n-2
        t = rpb*math.sqrt(df/((1.0-rpb+TINY)*(1.0+rpb+TINY)))
        prob = betai(0.5*df,0.5,df/(df+t*t))  # t already a float
        return rpb, prob


def lkendalltau(x,y):
    """
Calculates Kendall's tau ... correlation of ordinal data.  Adapted
from function kendl1 in Numerical Recipies.  Needs good test-routine.@@@

Usage:   lkendalltau(x,y)
Returns: Kendall's tau, two-tailed p-value
"""
    n1 = 0
    n2 = 0
    iss = 0
    for j in range(len(x)-1):
        for k in range(j,len(y)):
            a1 = x[j] - x[k]
            a2 = y[j] - y[k]
            aa = a1 * a2
            if (aa):             # neither list has a tie
                n1 = n1 + 1
                n2 = n2 + 1
                if aa > 0:
                    iss = iss + 1
                else:
                    iss = iss -1
            else:
                if (a1):
                    n1 = n1 + 1
                else:
                    n2 = n2 + 1
    tau = iss / math.sqrt(n1*n2)
    svar = (4.0*len(x)+10.0) / (9.0*len(x)*(len(x)-1))
    z = tau / math.sqrt(svar)
    prob = erfcc(abs(z)/1.4142136)
    return tau, prob


def llinregress(x,y):
    """
Calculates a regression line on x,y pairs.

Usage:   llinregress(x,y)     x,y are equal-length lists of x-y coordinates
Returns: slope, intercept, r, two-tailed prob, sterr-of-estimate
"""
    TINY = 1.0e-20
    if len(x) <> len(y):
        raise ValueError, 'Input values not paired in linregress.  Aborting.'
    n = len(x)
    x = map(float,x)
    y = map(float,y)
    xmean = mean(x)
    ymean = mean(y)
    r_num = float(n*(summult(x,y)) - sum(x)*sum(y))
    r_den = math.sqrt((n*ss(x) - square_of_sums(x))*(n*ss(y)-square_of_sums(y)))
    r = r_num / r_den
    z = 0.5*math.log((1.0+r+TINY)/(1.0-r+TINY))
    df = n-2
    t = r*math.sqrt(df/((1.0-r+TINY)*(1.0+r+TINY)))
    prob = betai(0.5*df,0.5,df/(df+t*t))
    slope = r_num / float(n*ss(x) - square_of_sums(x))
    intercept = ymean - slope*xmean
    sterrest = math.sqrt(1-r*r)*samplestdev(y)
    return slope, intercept, r, prob, sterrest


####################################
#####  INFERENTIAL STATISTICS  #####
####################################

def lttest_1samp(a,popmean,printit=0,name='Sample',writemode='a'):
    """
Calculates the t-obtained for the independent samples T-test on ONE group
of scores a, given a population mean.  If printit=1, results are printed
to the screen.  If printit='filename', the results are output to 'filename'
using the given writemode (default=append).  Returns t-value, and prob.

Usage:   lttest_1samp(a,popmean,Name='Sample',printit=0,writemode='a')
Returns: t-value, two-tailed prob
"""
    x = mean(a)
    v = var(a)
    n = len(a)
    df = n-1
    svar = ((n-1)*v)/float(df)
    t = (x-popmean)/math.sqrt(svar*(1.0/n))
    prob = betai(0.5*df,0.5,float(df)/(df+t*t))

    if printit <> 0:
        statname = 'Single-sample T-test.'
        outputpairedstats(printit,writemode,
                          'Population','--',popmean,0,0,0,
                          name,n,x,v,min(a),max(a),
                          statname,t,prob)
    return t,prob


def lttest_ind (a, b, printit=0, name1='Samp1', name2='Samp2', writemode='a'):
    """
Calculates the t-obtained T-test on TWO INDEPENDENT samples of
scores a, and b.  From Numerical Recipies, p.483.  If printit=1, results
are printed to the screen.  If printit='filename', the results are output
to 'filename' using the given writemode (default=append).  Returns t-value,
and prob.

Usage:   lttest_ind(a,b,printit=0,name1='Samp1',name2='Samp2',writemode='a')
Returns: t-value, two-tailed prob
"""
    x1 = mean(a)
    x2 = mean(b)
    v1 = stdev(a)**2
    v2 = stdev(b)**2
    n1 = len(a)
    n2 = len(b)
    df = n1+n2-2
    svar = ((n1-1)*v1+(n2-1)*v2)/float(df)
    t = (x1-x2)/math.sqrt(svar*(1.0/n1 + 1.0/n2))
    prob = betai(0.5*df,0.5,df/(df+t*t))

    if printit <> 0:
        statname = 'Independent samples T-test.'
        outputpairedstats(printit,writemode,
                          name1,n1,x1,v1,min(a),max(a),
                          name2,n2,x2,v2,min(b),max(b),
                          statname,t,prob)
    return t,prob


def lttest_rel (a,b,printit=0,name1='Sample1',name2='Sample2',writemode='a'):
    """
Calculates the t-obtained T-test on TWO RELATED samples of scores,
a and b.  From Numerical Recipies, p.483.  If printit=1, results are
printed to the screen.  If printit='filename', the results are output to
'filename' using the given writemode (default=append).  Returns t-value,
and prob.

Usage:   lttest_rel(a,b,printit=0,name1='Sample1',name2='Sample2',writemode='a')
Returns: t-value, two-tailed prob
"""
    if len(a)<>len(b):
        raise ValueError, 'Unequal length lists in ttest_rel.'
    x1 = mean(a)
    x2 = mean(b)
    v1 = var(a)
    v2 = var(b)
    n = len(a)
    cov = 0
    for i in range(len(a)):
        cov = cov + (a[i]-x1) * (b[i]-x2)
    df = n-1
    cov = cov / float(df)
    sd = math.sqrt((v1+v2 - 2.0*cov)/float(n))
    t = (x1-x2)/sd
    prob = betai(0.5*df,0.5,df/(df+t*t))

    if printit <> 0:
        statname = 'Related samples T-test.'
        outputpairedstats(printit,writemode,
                          name1,n,x1,v1,min(a),max(a),
                          name2,n,x2,v2,min(b),max(b),
                          statname,t,prob)
    return t, prob


def lchisquare(f_obs,f_exp=None):
    """
Calculates a one-way chi square for list of observed frequencies and returns
the result.  If no expected frequencies are given, the total N is assumed to
be equally distributed across all groups.

Usage:   lchisquare(f_obs, f_exp=None)   f_obs = list of observed cell freq.
Returns: chisquare-statistic, associated p-value
"""
    k = len(f_obs)               # number of groups
    if f_exp == None:
        f_exp = [sum(f_obs)/float(k)] * len(f_obs) # create k bins with = freq.
    chisq = 0
    for i in range(len(f_obs)):
        chisq = chisq + (f_obs[i]-f_exp[i])**2 / float(f_exp[i])
    return chisq, chisqprob(chisq, k-1)


def lks_2samp (data1,data2):
    """
Computes the Kolmogorov-Smirnof statistic on 2 samples.  From
Numerical Recipies in C, page 493.

Usage:   lks_2samp(data1,data2)   data1&2 are lists of values for 2 conditions
Returns: KS D-value, associated p-value
"""
    j1 = 0
    j2 = 0
    fn1 = 0.0
    fn2 = 0.0
    n1 = len(data1)
    n2 = len(data2)
    en1 = n1
    en2 = n2
    d = 0.0
    data1.sort()
    data2.sort()
    while j1 < n1 and j2 < n2:
        d1=data1[j1]
        d2=data2[j2]
        if d1 <= d2:
            fn1 = (j1)/float(en1)
            j1 = j1 + 1
        if d2 <= d1:
            fn2 = (j2)/float(en2)
            j2 = j2 + 1
        dt = (fn2-fn1)
        if math.fabs(dt) > math.fabs(d):
            d = dt
    try:
        en = math.sqrt(en1*en2/float(en1+en2))
        prob = ksprob((en+0.12+0.11/en)*abs(d))
    except:
        prob = 1.0
    return d, prob


def lmannwhitneyu(x,y):
    """
Calculates a Mann-Whitney U statistic on the provided scores and
returns the result.  Use only when the n in each condition is < 20 and
you have 2 independent samples of ranks.  NOTE: Mann-Whitney U is
significant if the u-obtained is LESS THAN or equal to the critical
value of U found in the tables.  Equivalent to Kruskal-Wallis H with
just 2 groups.

Usage:   lmannwhitneyu(data)
Returns: u-statistic, one-tailed p-value (i.e., p(z(U)))
"""
    n1 = len(x)
    n2 = len(y)
    ranked = rankdata(x+y)
    rankx = ranked[0:n1]       # get the x-ranks
    ranky = ranked[n1:]     # the rest are y-ranks
    u1 = n1*n2 + (n1*(n1+1))/2.0 - sum(rankx)  # calc U for x
    u2 = n1*n2 - u1                         # remainder is U for y
    bigu = max(u1,u2)
    smallu = min(u1,u2)
    T = math.sqrt(tiecorrect(ranked))  # correction factor for tied scores
    if T == 0:
        raise ValueError, 'All numbers are identical in lmannwhitneyu'
    sd = math.sqrt(T*n1*n2*(n1+n2+1)/12.0)
    z = abs((bigu-n1*n2/2.0) / sd)  # normal approximation for prob calc
    return smallu, 1.0 - zprob(z)


def ltiecorrect(rankvals):
    """
Corrects for ties in Mann Whitney U and Kruskal Wallis H tests.  See
Siegel, S. (1956) Nonparametric Statistics for the Behavioral Sciences.
New York: McGraw-Hill.  Code adapted from |Stat rankind.c code.

Usage:   ltiecorrect(rankvals)
Returns: T correction factor for U or H
"""
    sorted,posn = shellsort(rankvals)
    n = len(sorted)
    T = 0.0
    i = 0
    while (i<n-1):
        if sorted[i] == sorted[i+1]:
            nties = 1
            while (i<n-1) and (sorted[i] == sorted[i+1]):
                nties = nties +1
                i = i +1
            T = T + nties**3 - nties
        i = i+1
    T = T / float(n**3-n)
    return 1.0 - T


def lranksums(x,y):
    """
Calculates the rank sums statistic on the provided scores and
returns the result.  Use only when the n in each condition is > 20 and you
have 2 independent samples of ranks.

Usage:   lranksums(x,y)
Returns: a z-statistic, two-tailed p-value
"""
    n1 = len(x)
    n2 = len(y)
    alldata = x+y
    ranked = rankdata(alldata)
    x = ranked[:n1]
    y = ranked[n1:]
    s = sum(x)
    expected = n1*(n1+n2+1) / 2.0
    z = (s - expected) / math.sqrt(n1*n2*(n1+n2+1)/12.0)
    prob = 2*(1.0 -zprob(abs(z)))
    return z, prob


def lwilcoxont(x,y):
    """
Calculates the Wilcoxon T-test for related samples and returns the
result.  A non-parametric T-test.

Usage:   lwilcoxont(x,y)
Returns: a t-statistic, two-tail probability estimate
"""
    if len(x) <> len(y):
        raise ValueError, 'Unequal N in wilcoxont.  Aborting.'
    d=[]
    for i in range(len(x)):
        diff = x[i] - y[i]
        if diff <> 0:
            d.append(diff)
    count = len(d)
    absd = map(abs,d)
    absranked = rankdata(absd)
    r_plus = 0.0
    r_minus = 0.0
    for i in range(len(absd)):
        if d[i] < 0:
            r_minus = r_minus + absranked[i]
        else:
            r_plus = r_plus + absranked[i]
    wt = min(r_plus, r_minus)
    mn = count * (count+1) * 0.25
    se =  math.sqrt(count*(count+1)*(2.0*count+1.0)/24.0)
    z = math.fabs(wt-mn) / se
    prob = 2*(1.0 -zprob(abs(z)))
    return wt, prob


def lkruskalwallish(*args):
    """
The Kruskal-Wallis H-test is a non-parametric ANOVA for 3 or more
groups, requiring at least 5 subjects in each group.  This function
calculates the Kruskal-Wallis H-test for 3 or more independent samples
and returns the result.

Usage:   lkruskalwallish(*args)
Returns: H-statistic (corrected for ties), associated p-value
"""
    args = list(args)
    n = [0]*len(args)
    all = []
    n = map(len,args)
    for i in range(len(args)):
        all = all + args[i]
    ranked = rankdata(all)
    T = tiecorrect(ranked)
    for i in range(len(args)):
        args[i] = ranked[0:n[i]]
        del ranked[0:n[i]]
    rsums = []
    for i in range(len(args)):
        rsums.append(sum(args[i])**2)
        rsums[i] = rsums[i] / float(n[i])
    ssbn = sum(rsums)
    totaln = sum(n)
    h = 12.0 / (totaln*(totaln+1)) * ssbn - 3*(totaln+1)
    df = len(args) - 1
    if T == 0:
        raise ValueError, 'All numbers are identical in lkruskalwallish'
    h = h / float(T)
    return h, chisqprob(h,df)


def lfriedmanchisquare(*args):
    """
Friedman Chi-Square is a non-parametric, one-way within-subjects
ANOVA.  This function calculates the Friedman Chi-square test for repeated
measures and returns the result, along with the associated probability
value.  It assumes 3 or more repeated measures.  Only 3 levels requires a
minimum of 10 subjects in the study.  Four levels requires 5 subjects per
level(??).

Usage:   lfriedmanchisquare(*args)
Returns: chi-square statistic, associated p-value
"""
    k = len(args)
    if k < 3:
        raise ValueError, 'Less than 3 levels.  Friedman test not appropriate.'
    n = len(args[0])
    data = apply(pstat.abut,tuple(args))
    for i in range(len(data)):
        data[i] = rankdata(data[i])
    ssbn = 0
    for i in range(k):
        ssbn = ssbn + sum(args[i])**2
    chisq = 12.0 / (k*n*(k+1)) * ssbn - 3*n*(k+1)
    return chisq, chisqprob(chisq,k-1)


####################################
####  PROBABILITY CALCULATIONS  ####
####################################

def lchisqprob(chisq,df):
    """
    Returns the (1-tailed) probability value associated with the provided
    chi-square value and df.  Adapted from chisq.c in Gary Perlman's |Stat.

    Usage:   lchisqprob(chisq,df)
    """
    BIG = 20.0
    def ex(x):
        BIG = 20.0
        if x < -BIG:
            return 0.0
        else:
            return math.exp(x)

    if chisq <=0 or df < 1:
        return 1.0
    a = 0.5 * chisq
    if df%2 == 0:
        even = 1
    else:
        even = 0
    if df > 1:
        y = ex(-a)
    if even:
        s = y
    else:
        s = 2.0 * zprob(-math.sqrt(chisq))
    if (df > 2):
        chisq = 0.5 * (df - 1.0)
        if even:
            z = 1.0
        else:
            z = 0.5
        if a > BIG:
            if even:
                e = 0.0
            else:
                e = math.log(math.sqrt(math.pi))
            c = math.log(a)
            while (z <= chisq):
                e = math.log(z) + e
                s = s + ex(c*z-a-e)
                z = z + 1.0
            return s
        else:
            if even:
                e = 1.0
            else:
                e = 1.0 / math.sqrt(math.pi) / math.sqrt(a)
            c = 0.0
            while (z <= chisq):
                e = e * (a/float(z))
                c = c + e
                z = z + 1.0
            return (c*y+s)
    else:
        return s


def lerfcc(x):
    """
    Returns the complementary error function erfc(x) with fractional
    error everywhere less than 1.2e-7.  Adapted from Numerical Recipies.

    Usage:   lerfcc(x)
    """
    z = abs(x)
    t = 1.0 / (1.0+0.5*z)
    ans = t * math.exp(-z*z-1.26551223 + t*(1.00002368+t*(0.37409196+t*(0.09678418+t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+t*(-0.82215223+t*0.17087277)))))))))
    if x >= 0:
        return ans
    else:
        return 2.0 - ans


def zprob(z):
    """
    Returns the area under the normal curve 'to the left of' the given z value.
    Thus,
        for z<0, zprob(z) = 1-tail probability
        for z>0, 1.0-zprob(z) = 1-tail probability
        for any z, 2.0*(1.0-zprob(abs(z))) = 2-tail probability
    Adapted from z.c in Gary Perlman's |Stat.

    Usage:   lzprob(z)
    """
    Z_MAX = 6.0 # maximum meaningful z-value
    if z == 0.0:
        x = 0.0
    else:
        y = 0.5 * math.fabs(z)
        if y >= (Z_MAX*0.5):
            x = 1.0
        elif (y < 1.0):
            w = y*y
            x = ((((((((0.000124818987 * w
                        -0.001075204047) * w +0.005198775019) * w
                      -0.019198292004) * w +0.059054035642) * w
                    -0.151968751364) * w +0.319152932694) * w
                  -0.531923007300) * w +0.797884560593) * y * 2.0
        else:
            y = y - 2.0
            x = (((((((((((((-0.000045255659 * y
                             +0.000152529290) * y -0.000019538132) * y
                           -0.000676904986) * y +0.001390604284) * y
                         -0.000794620820) * y -0.002034254874) * y
                       +0.006549791214) * y -0.010557625006) * y
                     +0.011630447319) * y -0.009279453341) * y
                   +0.005353579108) * y -0.002141268741) * y
                 +0.000535310849) * y +0.999936657524
    if z > 0.0:
        prob = ((x+1.0)*0.5)
    else:
        prob = ((1.0-x)*0.5)
    return prob


def lksprob(alam):
    """
    Computes a Kolmolgorov-Smirnov t-test significance level.  Adapted from
    Numerical Recipies.

    Usage:   lksprob(alam)
    """
    fac = 2.0
    sum = 0.0
    termbf = 0.0
    a2 = -2.0*alam*alam
    for j in range(1,201):
        term = fac*math.exp(a2*j*j)
        sum = sum + term
        if math.fabs(term) <= (0.001*termbf) or math.fabs(term) < (1.0e-8*sum):
            return sum
        fac = -fac
        termbf = math.fabs(term)
    return 1.0           # Get here only if fails to converge; was 0.0!!


def lfprob (dfnum, dfden, F):
    """
    Returns the (1-tailed) significance level (p-value) of an F
    statistic given the degrees of freedom for the numerator (dfR-dfF) and
    the degrees of freedom for the denominator (dfF).

    Usage:   lfprob(dfnum, dfden, F)   where usually dfnum=dfbn, dfden=dfwn
    """
    p = betai(0.5*dfden, 0.5*dfnum, dfden/float(dfden+dfnum*F))
    return p


def lbetacf(a,b,x):
    """
    This function evaluates the continued fraction form of the incomplete
    Beta function, betai.  (Adapted from: Numerical Recipies in C.)

    Usage:   lbetacf(a,b,x)
    """
    ITMAX = 200
    EPS = 3.0e-7

    bm = az = am = 1.0
    qab = a+b
    qap = a+1.0
    qam = a-1.0
    bz = 1.0-qab*x/qap
    for i in range(ITMAX+1):
        em = float(i+1)
        tem = em + em
        d = em*(b-em)*x/((qam+tem)*(a+tem))
        ap = az + d*am
        bp = bz+d*bm
        d = -(a+em)*(qab+em)*x/((qap+tem)*(a+tem))
        app = ap+d*az
        bpp = bp+d*bz
        aold = az
        am = ap/bpp
        bm = bp/bpp
        az = app/bpp
        bz = 1.0
        if (abs(az-aold)<(EPS*abs(az))):
            return az
    print 'a or b too big, or ITMAX too small in Betacf.'


def lgammln(xx):
    """
Returns the gamma function of xx.
    Gamma(z) = Integral(0,infinity) of t^(z-1)exp(-t) dt.
(Adapted from: Numerical Recipies in C.)

Usage:   lgammln(xx)
"""

    coeff = [76.18009173, -86.50532033, 24.01409822, -1.231739516,
             0.120858003e-2, -0.536382e-5]
    x = xx - 1.0
    tmp = x + 5.5
    tmp = tmp - (x+0.5)*math.log(tmp)
    ser = 1.0
    for j in range(len(coeff)):
        x = x + 1
        ser = ser + coeff[j]/x
    return -tmp + math.log(2.50662827465*ser)


def lbetai(a,b,x):
    """
Returns the incomplete beta function:

    I-sub-x(a,b) = 1/B(a,b)*(Integral(0,x) of t^(a-1)(1-t)^(b-1) dt)

where a,b>0 and B(a,b) = G(a)*G(b)/(G(a+b)) where G(a) is the gamma
function of a.  The continued fraction formulation is implemented here,
using the betacf function.  (Adapted from: Numerical Recipies in C.)

Usage:   lbetai(a,b,x)
"""
    if (x<0.0 or x>1.0):
        raise ValueError, 'Bad x in lbetai'
    if (x==0.0 or x==1.0):
        bt = 0.0
    else:
        bt = math.exp(gammln(a+b)-gammln(a)-gammln(b)+a*math.log(x)+b*
                      math.log(1.0-x))
    if (x<(a+1.0)/(a+b+2.0)):
        return bt*betacf(a,b,x)/float(a)
    else:
        return 1.0-bt*betacf(b,a,1.0-x)/float(b)


####################################
#######  ANOVA CALCULATIONS  #######
####################################

def lF_oneway(*lists):
    """
Performs a 1-way ANOVA, returning an F-value and probability given
any number of groups.  From Heiman, pp.394-7.

Usage:   F_oneway(*lists)   where *lists is any number of lists, one per
                                  treatment group
Returns: F value, one-tailed p-value
"""
    a = len(lists)         # ANOVA on 'a' groups, each in it's own list
    means = [0]*a
    vars = [0]*a
    ns = [0]*a
    alldata = []
    tmp = map(N.array,lists)
    means = map(amean,tmp)
    vars = map(avar,tmp)
    ns = map(len,lists)
    for i in range(len(lists)):
        alldata = alldata + lists[i]
    alldata = N.array(alldata)
    bign = len(alldata)
    sstot = ass(alldata)-(asquare_of_sums(alldata)/float(bign))
    ssbn = 0
    for list in lists:
        ssbn = ssbn + asquare_of_sums(N.array(list))/float(len(list))
    ssbn = ssbn - (asquare_of_sums(alldata)/float(bign))
    sswn = sstot-ssbn
    dfbn = a-1
    dfwn = bign - a
    msb = ssbn/float(dfbn)
    msw = sswn/float(dfwn)
    f = msb/msw
    prob = fprob(dfbn,dfwn,f)
    return f, prob


def lF_value (ER,EF,dfnum,dfden):
    """
Returns an F-statistic given the following:
        ER  = error associated with the null hypothesis (the Restricted model)
        EF  = error associated with the alternate hypothesis (the Full model)
        dfR-dfF = degrees of freedom of the numerator
        dfF = degrees of freedom associated with the denominator/Full model

Usage:   lF_value(ER,EF,dfnum,dfden)
"""
    return ((ER-EF)/float(dfnum) / (EF/float(dfden)))


####################################
########  SUPPORT FUNCTIONS  #######
####################################

def writecc (listoflists,file,writetype='w',extra=2):
    """
Writes a list of lists to a file in columns, customized by the max
size of items within the columns (max size of items in col, +2 characters)
to specified file.  File-overwrite is the default.

Usage:   writecc (listoflists,file,writetype='w',extra=2)
Returns: None
"""
    if type(listoflists[0]) not in [ListType,TupleType]:
        listoflists = [listoflists]
    outfile = open(file,writetype)
    rowstokill = []
    list2print = copy.deepcopy(listoflists)
    for i in range(len(listoflists)):
        if listoflists[i] == ['\n'] or listoflists[i]=='\n' or listoflists[i]=='dashes':
            rowstokill = rowstokill + [i]
    rowstokill.reverse()
    for row in rowstokill:
        del list2print[row]
    maxsize = [0]*len(list2print[0])
    for col in range(len(list2print[0])):
        items = pstat.colex(list2print,col)
        items = map(pstat.makestr,items)
        maxsize[col] = max(map(len,items)) + extra
    for row in listoflists:
        if row == ['\n'] or row == '\n':
            outfile.write('\n')
        elif row == ['dashes'] or row == 'dashes':
            dashes = [0]*len(maxsize)
            for j in range(len(maxsize)):
                dashes[j] = '-'*(maxsize[j]-2)
            outfile.write(pstat.lineincustcols(dashes,maxsize))
        else:
            outfile.write(pstat.lineincustcols(row,maxsize))
        outfile.write('\n')
    outfile.close()
    return None


def lincr(l,cap):       # to increment a list up to a max-list of 'cap'
    """
Simulate a counting system from an n-dimensional list.

Usage:   lincr(l,cap)   l=list to increment, cap=max values for each list pos'n
Returns: next set of values for list l, OR -1 (if overflow)
"""
    l[0] = l[0] + 1  # e.g., [0,0,0] --> [2,4,3] (=cap)
    for i in range(len(l)):
        if l[i] > cap[i] and i < len(l)-1: # if carryover AND not done
            l[i] = 0
            l[i+1] = l[i+1] + 1
        elif l[i] > cap[i] and i == len(l)-1: # overflow past last column, must be finished
            l = -1
    return l


def lsum (inlist):
    """
Returns the sum of the items in the passed list.

Usage:   lsum(inlist)
"""
    s = 0
    for item in inlist:
        s = s + item
    return s


def lcumsum (inlist):
    """
Returns a list consisting of the cumulative sum of the items in the
passed list.

Usage:   lcumsum(inlist)
"""
    newlist = copy.deepcopy(inlist)
    for i in range(1,len(newlist)):
        newlist[i] = newlist[i] + newlist[i-1]
    return newlist


def lss(inlist):
    """
Squares each value in the passed list, adds up these squares and
returns the result.

Usage:   lss(inlist)
"""
    ss = 0
    for item in inlist:
        ss = ss + item*item
    return ss


def lsummult (list1,list2):
    """
Multiplies elements in list1 and list2, element by element, and
returns the sum of all resulting multiplications.  Must provide equal
length lists.

Usage:   lsummult(list1,list2)
"""
    if len(list1) <> len(list2):
        raise ValueError, "Lists not equal length in summult."
    s = 0
    for item1,item2 in pstat.abut(list1,list2):
        s = s + item1*item2
    return s


def lsumdiffsquared(x,y):
    """
Takes pairwise differences of the values in lists x and y, squares
these differences, and returns the sum of these squares.

Usage:   lsumdiffsquared(x,y)
Returns: sum[(x[i]-y[i])**2]
"""
    sds = 0
    for i in range(len(x)):
        sds = sds + (x[i]-y[i])**2
    return sds


def lsquare_of_sums(inlist):
    """
Adds the values in the passed list, squares the sum, and returns
the result.

Usage:   lsquare_of_sums(inlist)
Returns: sum(inlist[i])**2
"""
    s = sum(inlist)
    return float(s)*s


def lshellsort(inlist):
    """
Shellsort algorithm.  Sorts a 1D-list.

Usage:   lshellsort(inlist)
Returns: sorted-inlist, sorting-index-vector (for original list)
"""
    n = len(inlist)
    svec = copy.deepcopy(inlist)
    ivec = range(n)
    gap = n/2   # integer division needed
    while gap >0:
        for i in range(gap,n):
            for j in range(i-gap,-1,-gap):
                while j>=0 and svec[j]>svec[j+gap]:
                    temp        = svec[j]
                    svec[j]  = svec[j+gap]
                    svec[j+gap] = temp
                    itemp      = ivec[j]
                    ivec[j]  = ivec[j+gap]
                    ivec[j+gap] = itemp
        gap = gap / 2  # integer division needed
# svec is now sorted inlist, and ivec has the order svec[i] = vec[ivec[i]]
    return svec, ivec


def lrankdata(inlist):
    """
Ranks the data in inlist, dealing with ties appropritely.  Assumes
a 1D inlist.  Adapted from Gary Perlman's |Stat ranksort.

Usage:   lrankdata(inlist)
Returns: a list of length equal to inlist, containing rank scores
"""
    n = len(inlist)
    svec, ivec = shellsort(inlist)
    sumranks = 0
    dupcount = 0
    newlist = [0]*n
    for i in range(n):
        sumranks = sumranks + i
        dupcount = dupcount + 1
        if i==n-1 or svec[i] <> svec[i+1]:
            averank = sumranks / float(dupcount) + 1
            for j in range(i-dupcount+1,i+1):
                newlist[ivec[j]] = averank
            sumranks = 0
            dupcount = 0
    return newlist


def outputpairedstats(fname,writemode,name1,n1,m1,se1,min1,max1,name2,n2,m2,se2,min2,max2,statname,stat,prob):
    """
Prints or write to a file stats for two groups, using the name, n,
mean, sterr, min and max for each group, as well as the statistic name,
its value, and the associated p-value.

Usage:   outputpairedstats(fname,writemode,
                           name1,n1,mean1,stderr1,min1,max1,
                           name2,n2,mean2,stderr2,min2,max2,
                           statname,stat,prob)
Returns: None
"""
    suffix = ''                    # for *s after the p-value
    try:
        x = prob.shape
        prob = prob[0]
    except:
        pass
    if  prob < 0.001:  suffix = '  ***'
    elif prob < 0.01:  suffix = '  **'
    elif prob < 0.05:  suffix = '  *'
    title = [['Name','N','Mean','SD','Min','Max']]
    lofl = title+[[name1,n1,round(m1,3),round(math.sqrt(se1),3),min1,max1],
                  [name2,n2,round(m2,3),round(math.sqrt(se2),3),min2,max2]]
    if type(fname)<>StringType or len(fname)==0:
        print
        print statname
        print
        pstat.printcc(lofl)
        print
        try:
            if stat.shape == ():
                stat = stat[0]
            if prob.shape == ():
                prob = prob[0]
        except:
            pass
        print 'Test statistic = ',round(stat,3),'   p = ',round(prob,3),suffix
        print
    else:
        file = open(fname,writemode)
        file.write('\n'+statname+'\n\n')
        file.close()
        writecc(lofl,fname,'a')
        file = open(fname,'a')
        try:
            if stat.shape == ():
                stat = stat[0]
            if prob.shape == ():
                prob = prob[0]
        except:
            pass
        file.write(pstat.list2string(['\nTest statistic = ',round(stat,4),'   p = ',round(prob,4),suffix,'\n\n']))
        file.close()
    return None


def lfindwithin (data):
    """
Returns an integer representing a binary vector, where 1=within-
subject factor, 0=between.  Input equals the entire data 2D list (i.e.,
column 0=random factor, column -1=measured values (those two are skipped).
Note: input data is in |Stat format ... a list of lists ("2D list") with
one row per measured value, first column=subject identifier, last column=
score, one in-between column per factor (these columns contain level
designations on each factor).  See also stats.anova.__doc__.

Usage:   lfindwithin(data)   data in |Stat format
"""

    numfact = len(data[0])-1
    withinvec = 0
    for col in range(1,numfact):
        examplelevel = pstat.unique(pstat.colex(data,col))[0]
        rows = pstat.linexand(data,col,examplelevel)  # get 1 level of this factor
        factsubjs = pstat.unique(pstat.colex(rows,0))
        allsubjs = pstat.unique(pstat.colex(data,0))
        if len(factsubjs) == len(allsubjs):  # fewer Ss than scores on this factor?
            withinvec = withinvec + (1 << col)
    return withinvec


#########################################################
#########################################################
####### DISPATCH LISTS AND TUPLES TO ABOVE FCNS #########
#########################################################
#########################################################

## CENTRAL TENDENCY:
geometricmean = Dispatch ( (lgeometricmean, (ListType, TupleType)), )
harmonicmean = Dispatch ( (lharmonicmean, (ListType, TupleType)), )
mean = Dispatch ( (lmean, (ListType, TupleType)), )
median = Dispatch ( (lmedian, (ListType, TupleType)), )
medianscore = Dispatch ( (lmedianscore, (ListType, TupleType)), )
mode = Dispatch ( (lmode, (ListType, TupleType)), )

## MOMENTS:
moment = Dispatch ( (lmoment, (ListType, TupleType)), )
variation = Dispatch ( (lvariation, (ListType, TupleType)), )
skew = Dispatch ( (lskew, (ListType, TupleType)), )
kurtosis = Dispatch ( (lkurtosis, (ListType, TupleType)), )
describe = Dispatch ( (ldescribe, (ListType, TupleType)), )

## FREQUENCY STATISTICS:
itemfreq = Dispatch ( (litemfreq, (ListType, TupleType)), )
scoreatpercentile = Dispatch ( (lscoreatpercentile, (ListType, TupleType)), )
percentileofscore = Dispatch ( (lpercentileofscore, (ListType, TupleType)), )
histogram = Dispatch ( (lhistogram, (ListType, TupleType)), )
cumfreq = Dispatch ( (lcumfreq, (ListType, TupleType)), )
relfreq = Dispatch ( (lrelfreq, (ListType, TupleType)), )

## VARIABILITY:
obrientransform = Dispatch ( (lobrientransform, (ListType, TupleType)), )
samplevar = Dispatch ( (lsamplevar, (ListType, TupleType)), )
samplestdev = Dispatch ( (lsamplestdev, (ListType, TupleType)), )
var = Dispatch ( (lvar, (ListType, TupleType)), )
stdev = Dispatch ( (lstdev, (ListType, TupleType)), )
sterr = Dispatch ( (lsterr, (ListType, TupleType)), )
sem = Dispatch ( (lsem, (ListType, TupleType)), )
z = Dispatch ( (lz, (ListType, TupleType)), )
zs = Dispatch ( (lzs, (ListType, TupleType)), )

## TRIMMING FCNS:
trimboth = Dispatch ( (ltrimboth, (ListType, TupleType)), )
trim1 = Dispatch ( (ltrim1, (ListType, TupleType)), )

## CORRELATION FCNS:
paired = Dispatch ( (lpaired, (ListType, TupleType)), )
pearsonr = Dispatch ( (lpearsonr, (ListType, TupleType)), )
spearmanr = Dispatch ( (lspearmanr, (ListType, TupleType)), )
pointbiserialr = Dispatch ( (lpointbiserialr, (ListType, TupleType)), )
kendalltau = Dispatch ( (lkendalltau, (ListType, TupleType)), )
linregress = Dispatch ( (llinregress, (ListType, TupleType)), )

## INFERENTIAL STATS:
ttest_1samp = Dispatch ( (lttest_1samp, (ListType, TupleType)), )
ttest_ind = Dispatch ( (lttest_ind, (ListType, TupleType)), )
ttest_rel = Dispatch ( (lttest_rel, (ListType, TupleType)), )
chisquare = Dispatch ( (lchisquare, (ListType, TupleType)), )
ks_2samp = Dispatch ( (lks_2samp, (ListType, TupleType)), )
mannwhitneyu = Dispatch ( (lmannwhitneyu, (ListType, TupleType)), )
ranksums = Dispatch ( (lranksums, (ListType, TupleType)), )
tiecorrect = Dispatch ( (ltiecorrect, (ListType, TupleType)), )
wilcoxont = Dispatch ( (lwilcoxont, (ListType, TupleType)), )
kruskalwallish = Dispatch ( (lkruskalwallish, (ListType, TupleType)), )
friedmanchisquare = Dispatch ( (lfriedmanchisquare, (ListType, TupleType)), )

## PROBABILITY CALCS:
chisqprob = Dispatch ( (lchisqprob, (IntType, FloatType)), )
zprob = Dispatch ( (lzprob, (IntType, FloatType)), )
ksprob = Dispatch ( (lksprob, (IntType, FloatType)), )
fprob = Dispatch ( (lfprob, (IntType, FloatType)), )
betacf = Dispatch ( (lbetacf, (IntType, FloatType)), )
betai = Dispatch ( (lbetai, (IntType, FloatType)), )
erfcc = Dispatch ( (lerfcc, (IntType, FloatType)), )
gammln = Dispatch ( (lgammln, (IntType, FloatType)), )

## ANOVA FUNCTIONS:
F_oneway = Dispatch ( (lF_oneway, (ListType, TupleType)), )
F_value = Dispatch ( (lF_value, (ListType, TupleType)), )

## SUPPORT FUNCTIONS:
incr = Dispatch ( (lincr, (ListType, TupleType)), )
sum = Dispatch ( (lsum, (ListType, TupleType)), )
cumsum = Dispatch ( (lcumsum, (ListType, TupleType)), )
ss = Dispatch ( (lss, (ListType, TupleType)), )
summult = Dispatch ( (lsummult, (ListType, TupleType)), )
square_of_sums = Dispatch ( (lsquare_of_sums, (ListType, TupleType)), )
sumdiffsquared = Dispatch ( (lsumdiffsquared, (ListType, TupleType)), )
shellsort = Dispatch ( (lshellsort, (ListType, TupleType)), )
rankdata = Dispatch ( (lrankdata, (ListType, TupleType)), )
findwithin = Dispatch ( (lfindwithin, (ListType, TupleType)), )

# sort and discard these dispatches:

## CENTRAL TENDENCY:
geometricmean = Dispatch ( (lgeometricmean, (ListType, TupleType)),
                        (ageometricmean, (N.ndarray,)) )
harmonicmean = Dispatch ( (lharmonicmean, (ListType, TupleType)),
                       (aharmonicmean, (N.ndarray,)) )
mean = Dispatch ( (lmean, (ListType, TupleType)),
               (amean, (N.ndarray,)) )
median = Dispatch ( (lmedian, (ListType, TupleType)),
                 (amedian, (N.ndarray,)) )
medianscore = Dispatch ( (lmedianscore, (ListType, TupleType)),
                      (amedianscore, (N.ndarray,)) )
mode = Dispatch ( (lmode, (ListType, TupleType)),
               (amode, (N.ndarray,)) )
tmean = Dispatch ( (atmean, (N.ndarray,)) )
tvar = Dispatch ( (atvar, (N.ndarray,)) )
tstdev = Dispatch ( (atstdev, (N.ndarray,)) )
tsem = Dispatch ( (atsem, (N.ndarray,)) )

## VARIATION:
moment = Dispatch ( (lmoment, (ListType, TupleType)),
                 (amoment, (N.ndarray,)) )
variation = Dispatch ( (lvariation, (ListType, TupleType)),
                    (avariation, (N.ndarray,)) )
skew = Dispatch ( (lskew, (ListType, TupleType)),
               (askew, (N.ndarray,)) )
kurtosis = Dispatch ( (lkurtosis, (ListType, TupleType)),
                   (akurtosis, (N.ndarray,)) )
describe = Dispatch ( (ldescribe, (ListType, TupleType)),
                   (adescribe, (N.ndarray,)) )

## DISTRIBUTION TESTS

skewtest = Dispatch ( (askewtest, (ListType, TupleType)),
                   (askewtest, (N.ndarray,)) )
kurtosistest = Dispatch ( (akurtosistest, (ListType, TupleType)),
                       (akurtosistest, (N.ndarray,)) )
normaltest = Dispatch ( (anormaltest, (ListType, TupleType)),
                     (anormaltest, (N.ndarray,)) )

## FREQUENCY STATS:
itemfreq = Dispatch ( (litemfreq, (ListType, TupleType)),
                   (aitemfreq, (N.ndarray,)) )
scoreatpercentile = Dispatch ( (lscoreatpercentile, (ListType, TupleType)),
                            (ascoreatpercentile, (N.ndarray,)) )
percentileofscore = Dispatch ( (lpercentileofscore, (ListType, TupleType)),
                             (apercentileofscore, (N.ndarray,)) )
histogram = Dispatch ( (lhistogram, (ListType, TupleType)),
                    (ahistogram, (N.ndarray,)) )
cumfreq = Dispatch ( (lcumfreq, (ListType, TupleType)),
                  (acumfreq, (N.ndarray,)) )
relfreq = Dispatch ( (lrelfreq, (ListType, TupleType)),
                  (arelfreq, (N.ndarray,)) )

## VARIABILITY:
obrientransform = Dispatch ( (lobrientransform, (ListType, TupleType)),
                          (aobrientransform, (N.ndarray,)) )
samplevar = Dispatch ( (lsamplevar, (ListType, TupleType)),
                    (asamplevar, (N.ndarray,)) )
samplestdev = Dispatch ( (lsamplestdev, (ListType, TupleType)),
                      (asamplestdev, (N.ndarray,)) )
signaltonoise = Dispatch( (asignaltonoise, (N.ndarray,)),)
var = Dispatch ( (lvar, (ListType, TupleType)),
              (avar, (N.ndarray,)) )
stdev = Dispatch ( (lstdev, (ListType, TupleType)),
                (astdev, (N.ndarray,)) )
sterr = Dispatch ( (lsterr, (ListType, TupleType)),
                (asterr, (N.ndarray,)) )
sem = Dispatch ( (lsem, (ListType, TupleType)),
              (asem, (N.ndarray,)) )
z = Dispatch ( (lz, (ListType, TupleType)),
            (az, (N.ndarray,)) )
zs = Dispatch ( (lzs, (ListType, TupleType)),
             (azs, (N.ndarray,)) )

## TRIMMING FCNS:
threshold = Dispatch( (athreshold, (N.ndarray,)),)
trimboth = Dispatch ( (ltrimboth, (ListType, TupleType)),
                   (atrimboth, (N.ndarray,)) )
trim1 = Dispatch ( (ltrim1, (ListType, TupleType)),
                (atrim1, (N.ndarray,)) )

## CORRELATION FCNS:
paired = Dispatch ( (lpaired, (ListType, TupleType)),
                 (apaired, (N.ndarray,)) )
lincc = Dispatch ( (llincc, (ListType, TupleType)),
                   (alincc, (N.ndarray,)) )
pearsonr = Dispatch ( (lpearsonr, (ListType, TupleType)),
                   (apearsonr, (N.ndarray,)) )
spearmanr = Dispatch ( (lspearmanr, (ListType, TupleType)),
                    (aspearmanr, (N.ndarray,)) )
pointbiserialr = Dispatch ( (lpointbiserialr, (ListType, TupleType)),
                         (apointbiserialr, (N.ndarray,)) )
kendalltau = Dispatch ( (lkendalltau, (ListType, TupleType)),
                     (akendalltau, (N.ndarray,)) )
linregress = Dispatch ( (llinregress, (ListType, TupleType)),
                     (alinregress, (N.ndarray,)) )

## INFERENTIAL STATS:
ttest_1samp = Dispatch ( (lttest_1samp, (ListType, TupleType)),
                      (attest_1samp, (N.ndarray,)) )
ttest_ind = Dispatch ( (lttest_ind, (ListType, TupleType)),
                    (attest_ind, (N.ndarray,)) )
ttest_rel = Dispatch ( (lttest_rel, (ListType, TupleType)),
                    (attest_rel, (N.ndarray,)) )
chisquare = Dispatch ( (lchisquare, (ListType, TupleType)),
                    (achisquare, (N.ndarray,)) )
ks_2samp = Dispatch ( (lks_2samp, (ListType, TupleType)),
                   (aks_2samp, (N.ndarray,)) )
mannwhitneyu = Dispatch ( (lmannwhitneyu, (ListType, TupleType)),
                       (amannwhitneyu, (N.ndarray,)) )
tiecorrect = Dispatch ( (ltiecorrect, (ListType, TupleType)),
                     (atiecorrect, (N.ndarray,)) )
ranksums = Dispatch ( (lranksums, (ListType, TupleType)),
                   (aranksums, (N.ndarray,)) )
wilcoxont = Dispatch ( (lwilcoxont, (ListType, TupleType)),
                    (awilcoxont, (N.ndarray,)) )
kruskalwallish = Dispatch ( (lkruskalwallish, (ListType, TupleType)),
                         (akruskalwallish, (N.ndarray,)) )
friedmanchisquare = Dispatch ( (lfriedmanchisquare, (ListType, TupleType)),
                            (afriedmanchisquare, (N.ndarray,)) )

## PROBABILITY CALCS:
chisqprob = Dispatch ( (achisqprob, (IntType, FloatType)),
                    (achisqprob, (N.ndarray,)) )
zprob = Dispatch ( (lzprob, (IntType, FloatType)),
                (azprob, (N.ndarray,)) )
ksprob = Dispatch ( (lksprob, (IntType, FloatType)),
                 (aksprob, (N.ndarray,)) )
fprob = Dispatch ( (lfprob, (IntType, FloatType)),
                (afprob, (N.ndarray,)) )
betacf = Dispatch ( (lbetacf, (IntType, FloatType)),
                 (abetacf, (N.ndarray,)) )
betai = Dispatch ( (lbetai, (IntType, FloatType)),
                (abetai, (N.ndarray,)) )
erfcc = Dispatch ( (lerfcc, (IntType, FloatType)),
                (aerfcc, (N.ndarray,)) )
gammln = Dispatch ( (lgammln, (IntType, FloatType)),
                 (agammln, (N.ndarray,)) )

## ANOVA FUNCTIONS:
F_oneway = Dispatch ( (lF_oneway, (ListType, TupleType)),
                   (aF_oneway, (N.ndarray,)) )
F_value = Dispatch ( (lF_value, (ListType, TupleType)),
                  (aF_value, (N.ndarray,)) )

## SUPPORT FUNCTIONS:
incr = Dispatch ( (lincr, (ListType, TupleType, N.ndarray)), )
sum = Dispatch ( (lsum, (ListType, TupleType)),
              (asum, (N.ndarray,)) )
cumsum = Dispatch ( (lcumsum, (ListType, TupleType)),
                 (acumsum, (N.ndarray,)) )
ss = Dispatch ( (lss, (ListType, TupleType)),
             (ass, (N.ndarray,)) )
summult = Dispatch ( (lsummult, (ListType, TupleType)),
                  (asummult, (N.ndarray,)) )
square_of_sums = Dispatch ( (lsquare_of_sums, (ListType, TupleType)),
                         (asquare_of_sums, (N.ndarray,)) )
sumdiffsquared = Dispatch ( (lsumdiffsquared, (ListType, TupleType)),
                         (asumdiffsquared, (N.ndarray,)) )
shellsort = Dispatch ( (lshellsort, (ListType, TupleType)),
                    (ashellsort, (N.ndarray,)) )
rankdata = Dispatch ( (lrankdata, (ListType, TupleType)),
                   (arankdata, (N.ndarray,)) )
findwithin = Dispatch ( (lfindwithin, (ListType, TupleType)),
                     (afindwithin, (N.ndarray,)) )

if __name__ == "__main__":
    help(gammln)
