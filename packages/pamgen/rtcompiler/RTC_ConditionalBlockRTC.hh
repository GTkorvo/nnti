#ifndef _CONDITIONALBLOCKRTC_H
#define _CONDITIONALBLOCKRTC_H

#include "RTC_BlockRTC.hh"
#include "RTC_LineRTC.hh"
#include "RTC_VariableRTC.hh"
#include "RTC_TokenizerRTC.hh"

#include <string>
#include <map>

namespace PG_RuntimeCompiler {

/** 
 * A ConditionalBlock represents a block of code that only executes if
 * certain conditions are met. An example would be an if block or else if
 * block. ConditionalBlock extends Block because it is a Block.
 */

class ConditionalBlock : public Block
{
 public:
  
  /**
   * Constructor -> The constructor contructs the parent Block, finds and
   *                creates the conditional statement, and tells parent to
   *                create its sub statements. 
   * 
   * @param vars  - A map of already active variables
   * @param lines - The array of strings that represent the lines of the code. 
   * @param errs  - A string containing the errors that have been generated by
   *                the compiling of lines. If errs is not empty, then the 
   *                program has not compiled succesfully
   */
  ConditionalBlock(std::map<std::string, Variable*> vars, Tokenizer& lines, 
		   std::string& errs);

  /**
   * Destructor -> The destructor deletes the conditional statement
   */
  ~ConditionalBlock();

  /**
   * wasExecuted -> This method returns whether this block's been executed yet
   */
  bool wasExecuted() const {return _executed;}

  /**
   * reset -> This method resets _executed to its default, false. 
   *          This is necessary for ConditionalBlocks that are within loops. 
   */
  void reset() {_executed = false;}

  /**
   * execute -> This method executes this ConditionalBlock. The condition 
   *            statement is evaluated and the block is either entered or not 
   *            depending on this evaluation.
   */
  Value* execute();

 private:

  bool _executed; /**!< A bool that tells us if the Block was executed or the 
                   *    condition failed.
                   */
  
  Line* _condition; /**!< The Line of code that must be evaluated to be true 
                     *    if the ConditionalBlock is to be executed.
                     */
};

}
#endif
