/*
 * Copyright 2001 Vrije Universiteit, The Netherlands.
 * For full copyright and restrictions on use see the file COPYRIGHT in the
 * top level of the CCJ distribution.
 */

package CCJ;

import java.io.Serializable;

public interface Partitionable /* extends Serializable */{

    public void setElementAt(int index, int groupSize, Serializable object);
    public Serializable elementAt(int index, int groupSize);

}
