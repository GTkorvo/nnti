/*
 * Copyright 2001 Vrije Universiteit, The Netherlands.
 * For full copyright and restrictions on use see the file COPYRIGHT in the
 * top level of the CCJ distribution.
 */

package CCJ;

import java.io.Serializable;

final class MyEntry {
    MyEntry	next;
    int		tag;
    int		source;
    Serializable data;
    Serializable reply;
    boolean	received;
    boolean	replied;
    boolean	confirm;
    boolean	in_use;
    String	clearer;
    
    public String toString() {
	return "<tag=" + tag + ",source=" + source + ">";
    }
}
