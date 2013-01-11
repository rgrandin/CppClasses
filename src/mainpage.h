/**
 * @file mainpage.h
 * @mainpage Robert's General-Use Classes
 * @author 	Robert Grandin
 * @date 15 October 2010
 *
 * @section Introduction
 *
 * This documentation is for a collection of general-use C++ classes that I have developed
 * for use within my applications.  Most classes are intended to assist with the handling
 * of multi-dimensional datasets, but other classes are also included to handle linked-lists
 * and provide numerical routines.
 *
 * Whenever possible, these classes have been templated to provide the greatest-possible
 * flexibility and reusability of the code.  These classes are intended for numeric datatypes,
 * however, and unexpected behavior may occur if used on non-numeric types.
 *
 * Classes which simply store data, such as the Array classes, should work correctly for both
 * integer and floating-point numeric types.  Please note that some member callbacks of such
 * functions may produce unexpected results when used on an integer-type object.  For example,
 * the Mean() function of ArrayBase returns a value of a type which matches the templated type
 * of the class (e.g., for an array of integers it will return an integer value).  Since the
 * mean value of an array of integers is in-general not an integer, an alternative function,
 * MeanFloat() is provided to allow the user to capture single-precision decimal accuracy.
 * Such floating-point alternatives are not present in all cases where they might be useful,
 * but can be added as the need arises.
 *
 * Classes which exist primarily to manipulate data, such as the ODE, RootFinding, and Stats classes,
 * are intended to be templated by floating-point variable types.  Using integer versions of such
 * classes may produce unexpected and erroneous behavior.
 *
 *
 *
 *
 * @section Revisions
 *
 * @date 15 October 2010
 *	- Creation date.
 *
 * Creation and revision dates for individual classes are noted in the documentation
 * for those classes.  The creation date of 15 October 2010 is for the first versions
 * of my initial set of classes.
 *
 *
 *
 *
 * @section License
 *
 * Copyright (c) 2010-2012, Robert Grandin
 * All rights reserved.
 *
 * Redistribution and use of this file is permitted provided that the following
 * conditions are met:
 * 	-# 	Redistributions must produce the above copyright notice, this list of
 * 		conditions, and the following disclaimer in the documentation and/or
 * 		other materials provided with the distribution.
 * 	-#	Neither the name of the organization nor the names of its contributors
 * 		may be used to endorse or promote products derived from this software
 * 		without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY COPYRIGHT HOLDER "AS IS" AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO
 * EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING BUT NOT
 * LIMITING TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *
 *
 */

