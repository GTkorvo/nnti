// @HEADER
// ***********************************************************************
//
//               Java Implementation of the Petra Library
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

package Jpetra;

import java.io.PrintStream;

/*
 * Class.java
 *
 * Created on September 8, 2003, 2:06 PM
 */

/**
 *
 * @author  Jason Cross
 */
public class Output {
    private String label;
    private boolean rootPrint;
    private boolean vnodesPrint;
    private PrintStream rootStream;
    private PrintStream vnodesStream;
    
    public Output(String label, boolean rootPrint, PrintStream rootStream, boolean vnodesPrint, PrintStream vnodesStream) {
        this.label = label;
        this.rootPrint = rootPrint;
        this.rootStream = rootStream;
        this.vnodesPrint = vnodesPrint;
        this.vnodesStream = vnodesStream;
    }
    
    public void rootPrint(String message) {
        if (this.rootPrint) {
            rootStream.print(this.label + message);
        }
    }
    
    public void vnodesPrint(String message) {
        if (this.vnodesPrint) {
            vnodesStream.print(this.label + message);
        }
    }
    
    public void rootPrintln(String message) {
        if (this.rootPrint) {
            rootStream.println(this.label + message);
        }
    }
    
    public void vnodesPrintln(String message) {
        if (this.vnodesPrint) {
            vnodesStream.println(this.label + message);
        }
    }
    
    public boolean isRootPrint() {
        return this.rootPrint;
    }
    
    public boolean isVnodesPrint() {
        return this.vnodesPrint;
    }
    
    public void setRootPrint(boolean rootPrint) {
        this.rootPrint = rootPrint;
    }
    
    public void setVnodesPrint(boolean vnodesPrint) {
        this.vnodesPrint = vnodesPrint;
    }
    
    public void setRootStream(PrintStream rootStream) {
        this.rootStream = rootStream;
    }
    
    public void setVnodesStream(PrintStream vnodesStream) {
        this.vnodesStream = vnodesStream;
    }
}
