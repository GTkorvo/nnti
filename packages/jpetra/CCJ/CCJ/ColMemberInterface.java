/*
 * Copyright 2001 Vrije Universiteit, The Netherlands.
 * For full copyright and restrictions on use see the file COPYRIGHT in the
 * top level of the CCJ distribution.
 */

package CCJ;

import java.io.Serializable;

public interface ColMemberInterface {

    public void begin() throws CCJException;

    public Serializable broadcast(ColGroup group,
				  Serializable objectToDistribute,
				  int root)
			throws CCJException,
			       NoSuchMemberException;

    public Serializable reduce(ColGroup group,
			       Serializable dataObject,
			       Reducible reductionObject,
			       int root)
			throws CCJException,
			       NoSuchMemberException;

    public Serializable allReduce(ColGroup group,
				  Serializable dataObject,
				  Reducible reductionObject)
			throws CCJException,
			       NoSuchMemberException;

    public Partitionable flatGather(ColGroup group,
				    Partitionable rootObject,
				    Serializable dataObject,
				    int root)
			throws CCJException,
			       NoSuchMemberException;

    public Serializable flatScatter(ColGroup group,
				    Partitionable rootObject,
				    int root)
			throws CCJException,
			       NoSuchMemberException;

    public Partitionable gather(ColGroup group,
				Partitionable rootObject,
				Serializable dataObject,
				int root)
			throws CCJException,
			       NoSuchMemberException;

    public Partitionable allGather(ColGroup group,
				   Partitionable rootObject,
				   Serializable dataObject)
			throws CCJException,
			       NoSuchMemberException;

    public Serializable scatter(ColGroup group,
				Partitionable rootObject,
				int root)
			throws CCJException,
			       NoSuchMemberException;

    public void barrier(ColGroup group)
			throws CCJException, NoSuchMemberException;

    public void waitForCompletion();

    public int getMyRank(ColGroup group) throws CCJException;
}
