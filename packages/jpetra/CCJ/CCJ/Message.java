/*
 * Copyright 2001 Vrije Universiteit, The Netherlands.
 * For full copyright and restrictions on use see the file COPYRIGHT in the
 * top level of the CCJ distribution.
 */

package CCJ;

import java.io.Serializable;

final class Message implements Serializable {
    int tag;
    int source;
    Serializable data;

    transient boolean confirm;

    Message() {
	this.tag = -1;
	this.source = -1;
	this.data = null;
    }


    Message(int tag, int source, Serializable data) {
	this.tag = tag;
	this.source = source;
	this.data = data;
    }

}
