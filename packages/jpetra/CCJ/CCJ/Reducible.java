/*
 * Copyright 2001 Vrije Universiteit, The Netherlands.
 * For full copyright and restrictions on use see the file COPYRIGHT in the
 * top level of the CCJ distribution.
 */

package CCJ;

import java.io.Serializable;

public interface Reducible {

    public Serializable reduce(Serializable obj1, Serializable obj2);

}
