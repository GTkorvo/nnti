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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

/*
 * JpetraObject.java
 *
 * Created on June 4, 2001, 10:48 AM
 */

package Jpetra;
import java.io.Serializable;
import java.io.PrintStream;
import java.util.HashMap;
/**
 * JpetraObject is the base class for all classes in the
 * Jpetra package.
 *
 * @author  Michael William Boldt
 * @author Jason Cross
 * @version 
 */
public class JpetraObject extends java.lang.Object implements Serializable, Cloneable {
    private static HashMap outputStreams = new HashMap(4);
    /**
     * Creates a new JpetraObject
     */
    public JpetraObject() {
        if (outputStreams.isEmpty()) {
            outputStreams.put("STD", new Output("", true, System.out, false, System.out));
            outputStreams.put("ERR", new Output("Error: ", true, System.out, false, System.out));
            outputStreams.put("DEBUG", new Output("Debug: ", false, System.out, false, System.out));
            outputStreams.put("VERBOSE", new Output("", false, System.out, false, System.out));
        }
    }
    
    public static void setRootPrint (String key, boolean rootPrint) {
        Output out = (Output) outputStreams.get(key);
        out.setRootPrint(rootPrint);
    }

    public static void setRootStream (String key, PrintStream rootStream) {
        Output out = (Output) outputStreams.get(key);
        out.setRootStream(rootStream);
    }    
     
    public static void setVnodesPrint (String key, boolean vnodesPrint) {
        Output out = (Output) outputStreams.get(key);
        out.setVnodesPrint(vnodesPrint);
    }    
    
    public static void setVnodesStream (String key, PrintStream vnodesStream) {
        Output out = (Output) outputStreams.get(key);
        out.setVnodesStream(vnodesStream);
    } 
    
    public static void print (String key, String message) {
        Output out = (Output) outputStreams.get(key);
        out.print(message);
    }
 
    public static void println (String key, String message) {
        Output out = (Output) outputStreams.get(key);
        out.println(message);
    }
    
}
