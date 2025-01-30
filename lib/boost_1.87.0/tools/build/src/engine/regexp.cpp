/*
 * regcomp and regexec -- regsub and regerror are elsewhere
 *
 *  Copyright (c) 1986 by University of Toronto.
 *  Written by Henry Spencer.  Not derived from licensed software.
 *
 *  Permission is granted to anyone to use this software for any
 *  purpose on any computer system, and to redistribute it freely,
 *  subject to the following restrictions:
 *
 *  1. The author is not responsible for the consequences of use of
 *      this software, no matter how awful, even if they arise
 *      from defects in it.
 *
 *  2. The origin of this software must not be misrepresented, either
 *      by explicit claim or by omission.
 *
 *  3. Altered versions must be plainly marked as such, and must not
 *      be misrepresented as being the original software.
 *
 *** THIS IS AN ALTERED VERSION.  It was altered by John Gilmore,
 *** hoptoad!gnu, on 27 Dec 1986, to add \n as an alternative to |
 *** to assist in implementing egrep.
 *** THIS IS AN ALTERED VERSION.  It was altered by John Gilmore,
 *** hoptoad!gnu, on 27 Dec 1986, to add \< and \> for word-matching
 *** as in BSD grep and ex.
 *** THIS IS AN ALTERED VERSION.  It was altered by John Gilmore,
 *** hoptoad!gnu, on 28 Dec 1986, to optimize characters quoted with \.
 *** THIS IS AN ALTERED VERSION.  It was altered by James A. Woods,
 *** ames!jaw, on 19 June 1987, to quash a regcomp() redundancy.
 *** THIS IS AN ALTERED VERSION.  It was altered by Christopher Seiwald
 *** seiwald@vix.com, on 28 August 1993, for use in jam.  Regmagic.h
 *** was moved into regexp.h, and the include of regexp.h now uses "'s
 *** to avoid conflicting with the system regexp.h.  Const, bless its
 *** soul, was removed so it can compile everywhere.  The declaration
 *** of strchr() was in conflict on AIX, so it was removed (as it is
 *** happily defined in string.h).
 *** THIS IS AN ALTERED VERSION.  It was altered by Christopher Seiwald
 *** seiwald@perforce.com, on 20 January 2000, to use function prototypes.
 *** THIS IS AN ALTERED VERSION.  It was altered by René Ferdinand Rivera Morell
 *** on 2023-01, to convert to C++ and to be thread safe.
 *
 * Beware that some of this code is subtly aware of the way operator precedence
 * is structured in regular expressions. Serious changes in regular-expression
 * syntax might require a total rethink.
 */

#include "regexp.h"
#include "jam.h"
#include "output.h"
#include "strview.h"

#include <ctype.h>
#include <stdio.h>
#ifndef ultrix
#include <stdlib.h>
#endif
#include <string.h>

#include <string>
#include <unordered_map>

namespace b2 { namespace regex {

/*
 * The first byte of the regexp internal "program" is actually this magic
 * number; the start node begins in the second byte.
 */
#define MAGIC 0234

void regerror(char const * s);

/*
 * The "internal use only" fields in regexp.h are present to pass info from
 * compile to execute that permits the execute phase to run lots faster on
 * simple cases. They are:
 *
 *  regstart char that must begin a match; '\0' if none obvious.
 *  reganch  is the match anchored (at beginning-of-line only)?
 *  regmust  string (pointer into program) that match must include, or NULL.
 *  regmlen  length of regmust string.
 *
 * Regstart and reganch permit very fast decisions on suitable starting points
 * for a match, cutting down the work a lot.  Regmust permits fast rejection of
 * lines that cannot possibly match.  The regmust tests are costly enough that
 * regcomp() supplies a regmust only if the r.e. contains something potentially
 * expensive (at present, the only such thing detected is * or + at the start of
 * the r.e., which can involve a lot of backup). Regmlen is supplied because the
 * test in regexec() needs it and regcomp() is computing it anyway.
 */

/*
 * Structure for regexp "program". This is essentially a linear encoding of a
 * nondeterministic finite-state machine (aka syntax charts or "railroad normal
 * form" in parsing technology). Each node is an opcode plus a "next" pointer,
 * possibly plus an operand. "Next" pointers of all nodes except BRANCH
 * implement concatenation; a "next" pointer with a BRANCH on both ends of it is
 * connecting two alternatives. [Here we have one of the subtle syntax
 * dependencies: an individual BRANCH, as opposed to a collection of them, is
 * never concatenated with anything because of operator precedence.] The operand
 * of some types of node is a literal string; for others, it is a node leading
 * into a sub-FSM. In particular, the operand of a BRANCH node is the first node
 * of the branch. [NB this is *not* a tree structure: the tail of the branch
 * connects to the thing following the set of BRANCHes.] The opcodes are:
 */

/* definition   number  opnd?   meaning */
#define END 0 /* no   End of program. */
#define BOL 1 /* no   Match "" at beginning of line. */
#define EOL 2 /* no   Match "" at end of line. */
#define ANY 3 /* no   Match any one character. */
#define ANYOF 4 /* str  Match any character in this string. */
#define ANYBUT 5 /* str  Match any character not in this string. */
#define BRANCH 6 /* node Match this alternative, or the next... */
#define BACK 7 /* no   Match "", "next" ptr points backward. */
#define EXACTLY 8 /* str  Match this string. */
#define NOTHING 9 /* no   Match empty string. */
#define STAR 10 /* node Match this (simple) thing 0 or more times. */
#define PLUS 11 /* node Match this (simple) thing 1 or more times. */
#define WORDA 12 /* no   Match "" at wordchar, where prev is nonword */
#define WORDZ 13 /* no   Match "" at nonwordchar, where prev is word */
#define OPEN 20 /* no   Mark this point in input as start of #n. */
/*  OPEN+1 is number 1, etc. */
#define CLOSE 30 /* no   Analogous to OPEN. */

/*
 * Opcode notes:
 *
 * BRANCH   The set of branches constituting a single choice are hooked
 *      together with their "next" pointers, since precedence prevents
 *      anything being concatenated to any individual branch.  The
 *      "next" pointer of the last BRANCH in a choice points to the
 *      thing following the whole choice.  This is also where the
 *      final "next" pointer of each individual branch points; each
 *      branch starts with the operand node of a BRANCH node.
 *
 * BACK     Normal "next" pointers all implicitly point forward; BACK
 *      exists to make loop structures possible.
 *
 * STAR,PLUS    '?', and complex '*' and '+', are implemented as circular
 *      BRANCH structures using BACK.  Simple cases (one character
 *      per match) are implemented with STAR and PLUS for speed
 *      and to minimize recursive plunges.
 *
 * OPEN,CLOSE   ...are numbered at compile time.
 */

/*
 * A node is one char of opcode followed by two chars of "next" pointer.
 * "Next" pointers are stored as two 8-bit pieces, high order first.  The
 * value is a positive offset from the opcode of the node containing it.
 * An operand, if any, simply follows the node.  (Note that much of the
 * code generation knows about this implicit relationship.)
 *
 * Using two bytes for the "next" pointer is vast overkill for most things,
 * but allows patterns to get big without disasters.
 */
#define OP(p) (*(p))
#define NEXT(p) (((*((p) + 1) & 0377) << 8) + (*((p) + 2) & 0377))
#define OPERAND(p) ((p) + 3)

/*
 * See regmagic.h for one further detail of program structure.
 */

/*
 * Utility definitions.
 */
#ifndef CHARBITS
#define UCHARAT(p) ((int32_t) * (const unsigned char *)(p))
#else
#define UCHARAT(p) ((int32_t) * (p)&CHARBITS)
#endif

#define FAIL(m) \
{ \
regerror(m); \
return (NULL); \
}
#define ISMULT(c) ((c) == '*' || (c) == '+' || (c) == '?')

/*
 * Flags to be passed up and down.
 */
#define HASWIDTH 01 /* Known never to match null string. */
#define SIMPLE 02 /* Simple enough to be STAR/PLUS operand. */
#define SPSTART 04 /* Starts with * or +. */
#define WORST 0 /* Worst case. */

namespace {
char regdummy = 0;
/*
 - regnext - dig the "next" pointer out of a node
 */
template <typename C>
inline C * regnext(C * p)
{
	int32_t offset;

	if (p == &regdummy) return (NULL);

	offset = NEXT(p);
	if (offset == 0) return (NULL);

	if (OP(p) == BACK)
		return (p - offset);
	else
		return (p + offset);
}
} // namespace

// The compiled regex program to match with.
struct regex_prog
{
	std::string regexpr; /* The not-compiled regex.  */
	char regstart = 0; /* Internal use only. */
	char reganch = 0; /* Internal use only. */
	const char * regmust = nullptr; /* Internal use only. */
	int32_t regmlen = 0; /* Internal use only. */
	std::size_t progsize = 0; // The size of the program.
	char program[1]; /* Unwarranted chumminess with compiler. */

	regex_prog(const regex_prog &) = delete;
	regex_prog(regex_prog &&) = delete;
	inline regex_prog() = default;
};

using regex_prog_ptr = ::b2::jam::unique_jptr<regex_prog>;

struct compiler
{
	/*
	 * Global work variables for regcomp().
	 */
	char * regparse; /* Input-scan pointer. */
	int32_t regnpar; /* () count. */
	char * regcode; /* Code-emit pointer; &regdummy = don't. */
	int32_t regsize; /* Code size. */

	/*
	 - regcomp - compile a regular expression into internal code
	 *
	 * We can't allocate space until we know how big the compiled form will be,
	 * but we can't compile it (and thus know how big it is) until we've got a
	 * place to put the code.  So we cheat:  we compile it twice, once with code
	 * generation turned off and size counting turned on, and once "for real".
	 * This also means that we don't allocate space until we are sure that the
	 * thing really will compile successfully, and we never have to move the
	 * code and thus invalidate pointers into it.  (Note that it has to be in
	 * one piece because free() must be able to free it all.)
	 *
	 * Beware that the optimization-preparation code in here knows about some
	 * of the structure of the compiled regexp.
	 */
	regex_prog * regcomp(const char * exp)
	{
		regex_prog * r;
		char * scan;
		char * longest;
		int32_t len;
		int32_t flags;

		if (exp == NULL) FAIL("NULL argument");

			/* First pass: determine size, legality. */
#ifdef notdef
		if (exp[0] == '.' && exp[1] == '*') exp += 2; /* aid grep */
#endif
		regparse = (char *)exp;
		regnpar = 1;
		regsize = 0;
		regcode = &regdummy;
		regc(MAGIC);
		if (reg(0, &flags) == NULL) return (NULL);

		/* Small enough for pointer-storage convention? */
		if (regsize >= 32767L) /* Probably could be 65535L. */
			FAIL("regexp too big");

		/* Allocate space. */
		r = (regex_prog *)BJAM_MALLOC(sizeof(regex_prog) + regsize);
		if (r == NULL) FAIL("out of space");
		b2::jam::ctor_ptr<regex_prog>(r);
		r->regexpr = exp;
		r->progsize = regsize;

		/* Second pass: emit code. */
		regparse = (char *)exp;
		regnpar = 1;
		regcode = r->program;
		regc(MAGIC);
		if (reg(0, &flags) == NULL) return (NULL);

		/* Dig out information for optimizations. */
		r->regstart = '\0'; /* Worst-case defaults. */
		r->reganch = 0;
		r->regmust = NULL;
		r->regmlen = 0;
		scan = r->program + 1; /* First BRANCH. */
		if (OP(regnext(scan)) == END)
		{ /* Only one top-level choice. */
			scan = OPERAND(scan);

			/* Starting-point info. */
			if (OP(scan) == EXACTLY)
				r->regstart = *OPERAND(scan);
			else if (OP(scan) == BOL)
				r->reganch++;

			/*
			 * If there's something expensive in the r.e., find the
			 * longest literal string that must appear and make it the
			 * regmust.  Resolve ties in favor of later strings, since
			 * the regstart check works with the beginning of the r.e.
			 * and avoiding duplication strengthens checking.  Not a
			 * strong reason, but sufficient in the absence of others.
			 */
			if (flags & SPSTART)
			{
				longest = NULL;
				len = 0;
				for (; scan != NULL; scan = regnext(scan))
					if (OP(scan) == EXACTLY
						&& static_cast<int32_t>(strlen(OPERAND(scan))) >= len)
					{
						longest = OPERAND(scan);
						len = static_cast<int32_t>(strlen(OPERAND(scan)));
					}
				r->regmust = longest;
				r->regmlen = len;
			}
		}

		return (r);
	}

	/*
	 - reg - regular expression, i.e. main body or parenthesized thing
	 *
	 * Caller must absorb opening parenthesis.
	 *
	 * Combining parenthesis handling with the base level of regular expression
	 * is a trifle forced, but the need to tie the tails of the branches to what
	 * follows makes it hard to avoid.
	 */
	char * reg(int32_t paren, /* Parenthesized? */
		int32_t * flagp)
	{
		char * ret;
		char * br;
		char * ender;
		int32_t parno = 0;
		int32_t flags;

		*flagp = HASWIDTH; /* Tentatively. */

		/* Make an OPEN node, if parenthesized. */
		if (paren)
		{
			if (regnpar >= NSUBEXP) FAIL("too many ()");
			parno = regnpar;
			regnpar++;
			ret = regnode(OPEN + parno);
		}
		else
			ret = NULL;

		/* Pick up the branches, linking them together. */
		br = regbranch(&flags);
		if (br == NULL) return (NULL);
		if (ret != NULL)
			regtail(ret, br); /* OPEN -> first. */
		else
			ret = br;
		if (!(flags & HASWIDTH)) *flagp &= ~HASWIDTH;
		*flagp |= flags & SPSTART;
		while (*regparse == '|' || *regparse == '\n')
		{
			regparse++;
			br = regbranch(&flags);
			if (br == NULL) return (NULL);
			regtail(ret, br); /* BRANCH -> BRANCH. */
			if (!(flags & HASWIDTH)) *flagp &= ~HASWIDTH;
			*flagp |= flags & SPSTART;
		}

		/* Make a closing node, and hook it on the end. */
		ender = regnode((paren) ? CLOSE + parno : END);
		regtail(ret, ender);

		/* Hook the tails of the branches to the closing node. */
		for (br = ret; br != NULL; br = regnext(br)) regoptail(br, ender);

		/* Check for proper termination. */
		if (paren && *regparse++ != ')')
		{
			FAIL("unmatched ()");
		}
		else if (!paren && *regparse != '\0')
		{
			if (*regparse == ')')
			{
				FAIL("unmatched ()");
			}
			else
				FAIL("junk on end"); /* "Can't happen". */
			/* NOTREACHED */
		}

		return (ret);
	}

	/*
	 - regbranch - one alternative of an | operator
	 *
	 * Implements the concatenation operator.
	 */
	char * regbranch(int32_t * flagp)
	{
		char * ret;
		char * chain;
		char * latest;
		int32_t flags;

		*flagp = WORST; /* Tentatively. */

		ret = regnode(BRANCH);
		chain = NULL;
		while (*regparse != '\0' && *regparse != ')' && *regparse != '\n'
			&& *regparse != '|')
		{
			latest = regpiece(&flags);
			if (latest == NULL) return (NULL);
			*flagp |= flags & HASWIDTH;
			if (chain == NULL) /* First piece. */
				*flagp |= flags & SPSTART;
			else
				regtail(chain, latest);
			chain = latest;
		}
		if (chain == NULL) /* Loop ran zero times. */
			(void)regnode(NOTHING);

		return (ret);
	}

	/*
	 - regpiece - something followed by possible [*+?]
	 *
	 * Note that the branching code sequences used for ? and the general cases
	 * of * and + are somewhat optimized:  they use the same NOTHING node as
	 * both the endmarker for their branch list and the body of the last branch.
	 * It might seem that this node could be dispensed with entirely, but the
	 * endmarker role is not redundant.
	 */
	char * regpiece(int32_t * flagp)
	{
		char * ret;
		char op;
		char * next;
		int32_t flags;

		ret = regatom(&flags);
		if (ret == NULL) return (NULL);

		op = *regparse;
		if (!ISMULT(op))
		{
			*flagp = flags;
			return (ret);
		}

		if (!(flags & HASWIDTH) && op != '?') FAIL("*+ operand could be empty");
		*flagp = (op != '+') ? (WORST | SPSTART) : (WORST | HASWIDTH);

		if (op == '*' && (flags & SIMPLE))
			reginsert(STAR, ret);
		else if (op == '*')
		{
			/* Emit x* as (x&|), where & means "self". */
			reginsert(BRANCH, ret); /* Either x */
			regoptail(ret, regnode(BACK)); /* and loop */
			regoptail(ret, ret); /* back */
			regtail(ret, regnode(BRANCH)); /* or */
			regtail(ret, regnode(NOTHING)); /* null. */
		}
		else if (op == '+' && (flags & SIMPLE))
			reginsert(PLUS, ret);
		else if (op == '+')
		{
			/* Emit x+ as x(&|), where & means "self". */
			next = regnode(BRANCH); /* Either */
			regtail(ret, next);
			regtail(regnode(BACK), ret); /* loop back */
			regtail(next, regnode(BRANCH)); /* or */
			regtail(ret, regnode(NOTHING)); /* null. */
		}
		else if (op == '?')
		{
			/* Emit x? as (x|) */
			reginsert(BRANCH, ret); /* Either x */
			regtail(ret, regnode(BRANCH)); /* or */
			next = regnode(NOTHING); /* null. */
			regtail(ret, next);
			regoptail(ret, next);
		}
		regparse++;
		if (ISMULT(*regparse)) FAIL("nested *?+");

		return (ret);
	}

	/*
	 - regatom - the lowest level
	 *
	 * Optimization:  gobbles an entire sequence of ordinary characters so that
	 * it can turn them into a single node, which is smaller to store and
	 * faster to run.  Backslashed characters are exceptions, each becoming a
	 * separate node; the code is simpler that way and it's not worth fixing.
	 */
	char * regatom(int32_t * flagp)
	{
		char * ret;
		int32_t flags;

		*flagp = WORST; /* Tentatively. */

		switch (*regparse++)
		{
			/* FIXME: these chars only have meaning at beg/end of pat? */
			case '^': ret = regnode(BOL); break;
			case '$': ret = regnode(EOL); break;
			case '.':
				ret = regnode(ANY);
				*flagp |= HASWIDTH | SIMPLE;
				break;
			case '[':
			{
				int32_t classr;
				int32_t classend;

				if (*regparse == '^')
				{ /* Complement of range. */
					ret = regnode(ANYBUT);
					regparse++;
				}
				else
					ret = regnode(ANYOF);
				if (*regparse == ']' || *regparse == '-') regc(*regparse++);
				while (*regparse != '\0' && *regparse != ']')
				{
					if (*regparse == '-')
					{
						regparse++;
						if (*regparse == ']' || *regparse == '\0')
							regc('-');
						else
						{
							classr = UCHARAT(regparse - 2) + 1;
							classend = UCHARAT(regparse);
							if (classr > classend + 1) FAIL("invalid [] range");
							for (; classr <= classend; classr++) regc(classr);
							regparse++;
						}
					}
					else
						regc(*regparse++);
				}
				regc('\0');
				if (*regparse != ']') FAIL("unmatched []");
				regparse++;
				*flagp |= HASWIDTH | SIMPLE;
			}
			break;
			case '(':
				ret = reg(1, &flags);
				if (ret == NULL) return (NULL);
				*flagp |= flags & (HASWIDTH | SPSTART);
				break;
			case '\0':
			case '|':
			case '\n':
			case ')':
				FAIL("internal urp"); /* Supposed to be caught earlier. */
				break;
			case '?':
			case '+':
			case '*': FAIL("?+* follows nothing"); break;
			case '\\':
				switch (*regparse++)
				{
					case '\0': FAIL("trailing \\"); break;
					case '<': ret = regnode(WORDA); break;
					case '>': ret = regnode(WORDZ); break;
					/* FIXME: Someday handle \1, \2, ... */
					default:
						/* Handle general quoted chars in exact-match routine */
						goto de_fault;
				}
				break;
			de_fault:
			default:
				/*
				 * Encode a string of characters to be matched exactly.
				 *
				 * This is a bit tricky due to quoted chars and due to
				 * '*', '+', and '?' taking the SINGLE char previous
				 * as their operand.
				 *
				 * On entry, the char at regparse[-1] is going to go
				 * into the string, no matter what it is.  (It could be
				 * following a \ if we are entered from the '\' case.)
				 *
				 * Basic idea is to pick up a good char in  ch  and
				 * examine the next char.  If it's *+? then we twiddle.
				 * If it's \ then we frozzle.  If it's other magic char
				 * we push  ch  and terminate the string.  If none of the
				 * above, we push  ch  on the string and go around again.
				 *
				 *  regprev  is used to remember where "the current char"
				 * starts in the string, if due to a *+? we need to back
				 * up and put the current char in a separate, 1-char, string.
				 * When  regprev  is NULL,  ch  is the only char in the
				 * string; this is used in *+? handling, and in setting
				 * flags |= SIMPLE at the end.
				 */
				{
					char * regprev;
					char ch;

					regparse--; /* Look at cur char */
					ret = regnode(EXACTLY);
					for (regprev = 0;;)
					{
						ch = *regparse++; /* Get current char */
						switch (*regparse)
						{ /* look at next one */

							default:
								regc(ch); /* Add cur to string */
								break;

							case '.':
							case '[':
							case '(':
							case ')':
							case '|':
							case '\n':
							case '$':
							case '^':
							case '\0':
							/* FIXME, $ and ^ should not always be magic */
							magic:
								regc(ch); /* dump cur char */
								goto done; /* and we are done */

							case '?':
							case '+':
							case '*':
								if (!regprev) /* If just ch in str, */
									goto magic; /* use it */
								/* End mult-char string one early */
								regparse = regprev; /* Back up parse */
								goto done;

							case '\\':
								regc(ch); /* Cur char OK */
								switch (regparse[1])
								{ /* Look after \ */
									case '\0':
									case '<':
									case '>':
										/* FIXME: Someday handle \1, \2, ... */
										goto done; /* Not quoted */
									default:
										/* Backup point is \, scan * point is
										 * after it. */
										regprev = regparse;
										regparse++;
										continue; /* NOT break; */
								}
						}
						regprev = regparse; /* Set backup point */
					}
				done:
					regc('\0');
					*flagp |= HASWIDTH;
					if (!regprev) /* One char? */
						*flagp |= SIMPLE;
				}
				break;
		}

		return (ret);
	}

	/*
	 - regnode - emit a node
	 */
	char * /* Location. */
		regnode(int32_t op)
	{
		char * ret;
		char * ptr;

		ret = regcode;
		if (ret == &regdummy)
		{
			regsize += 3;
			return (ret);
		}

		ptr = ret;
		*ptr++ = op;
		*ptr++ = '\0'; /* Null "next" pointer. */
		*ptr++ = '\0';
		regcode = ptr;

		return (ret);
	}

	/*
	 - regc - emit (if appropriate) a byte of code
	 */
	void regc(int32_t b)
	{
		if (regcode != &regdummy)
			*regcode++ = b;
		else
			regsize++;
	}

	/*
	 - reginsert - insert an operator in front of already-emitted operand
	 *
	 * Means relocating the operand.
	 */
	void reginsert(char op, char * opnd)
	{
		char * src;
		char * dst;
		char * place;

		if (regcode == &regdummy)
		{
			regsize += 3;
			return;
		}

		src = regcode;
		regcode += 3;
		dst = regcode;
		while (src > opnd) *--dst = *--src;

		place = opnd; /* Op node, where operand used to be. */
		*place++ = op;
		*place++ = '\0';
		*place++ = '\0';
	}

	/*
	 - regtail - set the next-pointer at the end of a node chain
	 */
	void regtail(char * p, char * val)
	{
		char * scan;
		char * temp;
		size_t offset;

		if (p == &regdummy) return;

		/* Find last node. */
		scan = p;
		for (;;)
		{
			temp = regnext(scan);
			if (temp == NULL) break;
			scan = temp;
		}

		if (OP(scan) == BACK)
			offset = scan - val;
		else
			offset = val - scan;
		*(scan + 1) = (offset >> 8) & 0377;
		*(scan + 2) = offset & 0377;
	}

	/*
	 - regoptail - regtail on operand of first argument; nop if operandless
	 */

	void regoptail(char * p, char * val)
	{
		/* "Operandless" and "op != BRANCH" are synonymous in practice. */
		if (p == NULL || p == &regdummy || OP(p) != BRANCH) return;
		regtail(OPERAND(p), val);
	}

}; // struct compiler

regex_prog_ptr regex_comp(const char * string)
{
	compiler c;
	return regex_prog_ptr(c.regcomp(string));
}

struct executor
{
	/*
	 * regexec and friends
	 */

	/*
	 * Global work variables for regexec().
	 */
	string_view reg_in; /* String-input. */
	const char * reg_bol; /* Beginning of input, for ^ check. */

	/*
	 - regexec - match a regexp against a string
	 */
	bool regexec(
		const regex_prog & prog, regex_expr & expr, const string_view & string)
	{
		string_view s;
		/* Check validity of program. */
		if (UCHARAT(prog.program) != MAGIC)
		{
			regerror("corrupted program");
			return false;
		}

		/* If there is a "must appear" string, look for it. */
		if (prog.regmust != NULL
			&& string.find(prog.regmust, 0, prog.regmlen) == string_view::npos)
			return false; /* Not present. */

		/* Mark beginning of line for ^ . */
		reg_bol = string.begin();

		/* Simplest case:  anchored match need be tried only once. */
		if (prog.reganch) return regtry(prog, expr, string);

		/* Messy cases:  unanchored match. */
		s = string;
		if (prog.regstart != '\0') /* We know what char it must start with. */
			for (auto j = s.find(prog.regstart); j != string_view::npos;
				 j = s.find(prog.regstart))
			{
				s = s.substr(j);
				if (regtry(prog, expr, s)) return true;
				s = s.substr(1);
			}
		else
			/* We do not -- general case. */
			for (; !s.empty(); s = s.substr(1))
			{
				if (regtry(prog, expr, s)) return true;
			}

		/* Failure. */
		return false;
	}

	/*
	 * regtry() - try match at specific point.
	 * success == true
	 */

	inline bool regtry(
		const regex_prog & prog, regex_expr & expr, const string_view & string)
	{
		reg_in = string;
		expr = regex_expr {};
		if (regmatch(prog.program + 1, expr))
		{
			expr.sub[0] = string_view(string.begin(), reg_in.begin());
			return true;
		}
		return false;
	}

	/*
	 * regmatch() - main matching routine.
	 *
	 * Conceptually the strategy is simple: check to see whether the current
	 * node matches, call self recursively to see whether the rest matches, and
	 * then act accordingly. In practice we make some effort to avoid recursion,
	 * in particular by going through "ordinary" nodes (that do not need to know
	 * whether the rest of the match failed) by a loop instead of by recursion.
	 */

	bool regmatch(const char * prog, regex_expr & expr)
	{
		const char * scan; /* Current node. */
		const char * next; /* Next node. */

		scan = prog;
		while (scan != NULL)
		{
			next = regnext(scan);

			switch (OP(scan))
			{
				case BOL:
					if (reg_in.begin() != reg_bol) return false;
					break;
				case EOL:
					if (!reg_in.empty()) return false;
					break;
				case WORDA:
					/* Must be looking at a letter, digit, or _ */
					if (!reg_in.empty() && (!isalnum(reg_in[0]))
						&& reg_in[0] != '_')
						return false;
					/* Prev must be BOL or nonword */
					if (reg_in.begin() > reg_bol
						&& (isalnum(reg_in.begin()[-1])
							|| reg_in.begin()[-1] == '_'))
						return false;
					break;
				case WORDZ:
					/* Must be looking at non letter, digit, or _ */
					if (reg_in.empty() || isalnum(reg_in[0])
						|| reg_in[0] == '_')
						return false;
					/* We don't care what the previous char was */
					break;
				case ANY:
					if (reg_in.empty()) return false;
					reg_in = reg_in.substr(1);
					break;
				case EXACTLY:
				{
					size_t len;
					const char * opnd;

					opnd = OPERAND(scan);
					/* Inline the first character, for speed. */
					if (reg_in.empty() || *opnd != reg_in[0]) return false;
					len = strlen(opnd);
					if (len > 1 && reg_in.compare(0, len, opnd) != 0)
						return false;
					reg_in = reg_in.substr(len);
				}
				break;
				case ANYOF:
					if (reg_in.empty()
						|| strchr(OPERAND(scan), reg_in[0]) == NULL)
						return false;
					reg_in = reg_in.substr(1);
					break;
				case ANYBUT:
					if (reg_in.empty()
						|| strchr(OPERAND(scan), reg_in[0]) != NULL)
						return false;
					reg_in = reg_in.substr(1);
					break;
				case NOTHING: break;
				case BACK: break;
				case OPEN + 1:
				case OPEN + 2:
				case OPEN + 3:
				case OPEN + 4:
				case OPEN + 5:
				case OPEN + 6:
				case OPEN + 7:
				case OPEN + 8:
				case OPEN + 9:
				{
					int32_t no;

					no = OP(scan) - OPEN;
					auto save = reg_in;

					if (regmatch(next, expr))
					{
						/*
						 * Don't set startp if some later
						 * invocation of the same parentheses
						 * already has.
						 */
						if (expr.sub[no].begin() != nullptr)
							expr.sub[no] = string_view(
								save.begin(), expr.sub[no].begin());
						return true;
					}
					else
						return false;
				}
				break;
				case CLOSE + 1:
				case CLOSE + 2:
				case CLOSE + 3:
				case CLOSE + 4:
				case CLOSE + 5:
				case CLOSE + 6:
				case CLOSE + 7:
				case CLOSE + 8:
				case CLOSE + 9:
				{
					int32_t no;

					no = OP(scan) - CLOSE;
					auto save = reg_in;

					if (regmatch(next, expr))
					{
						/*
						 * Don't set endp if some later
						 * invocation of the same parentheses
						 * already has.
						 */
						if (expr.sub[no].begin() == nullptr)
							expr.sub[no] = string_view(save.begin(), 0);
						return true;
					}
					else
						return false;
				}
				break;
				case BRANCH:
				{
					if (OP(next) != BRANCH) /* No choice. */
						next = OPERAND(scan); /* Avoid recursion. */
					else
					{
						do
						{
							auto save = reg_in;
							if (regmatch(OPERAND(scan), expr)) return true;
							reg_in = save;
							scan = regnext(scan);
						}
						while (scan != NULL && OP(scan) == BRANCH);
						return false;
						/* NOTREACHED */
					}
				}
				break;
				case STAR:
				case PLUS:
				{
					char nextch;
					int32_t no;
					int32_t min;

					/*
					 * Lookahead to avoid useless match attempts
					 * when we know what character comes next.
					 */
					nextch = '\0';
					if (OP(next) == EXACTLY) nextch = *OPERAND(next);
					min = (OP(scan) == STAR) ? 0 : 1;
					auto save = reg_in;
					no = regrepeat(OPERAND(scan));
					while (no >= min)
					{
						/* If it could work, try it. */
						if (nextch == '\0' ||
							(reg_in.empty() ? '\0' : reg_in[0]) == nextch)
							if (regmatch(next, expr)) return true;
						/* Couldn't or didn't -- back up. */
						no--;
						reg_in = save.substr((std::max)(no, 0));
					}
					return false;
				}
				break;
				case END:
					return true; /* Success! */
					break;
				default:
					regerror("memory corruption");
					return false;
					break;
			}

			scan = next;
		}

		/*
		 * We get here only if there's trouble -- normally "case END" is
		 * the terminating point.
		 */
		regerror("corrupted pointers");
		return false;
	}

	/*
	 - regrepeat - repeatedly match something simple, report how many
	 */
	int32_t regrepeat(const char * p)
	{
		int32_t count = 0;
		const char * opnd;

		auto scan = reg_in;
		opnd = OPERAND(p);
		switch (OP(p))
		{
			case ANY:
				count = int32_t(scan.length());
				scan = scan.substr(count);
				break;
			case EXACTLY:
				while (!scan.empty() && *opnd == scan[0])
				{
					count++;
					scan = scan.substr(1);
				}
				break;
			case ANYOF:
				while (!scan.empty() && strchr(opnd, scan[0]) != NULL)
				{
					count++;
					scan = scan.substr(1);
				}
				break;
			case ANYBUT:
				while (!scan.empty() && strchr(opnd, scan[0]) == NULL)
				{
					count++;
					scan = scan.substr(1);
				}
				break;
			default: /* Oh dear.  Called inappropriately. */
				regerror("internal foulup");
				count = 0; /* Best compromise. */
				break;
		}
		reg_in = scan;

		return (count);
	}

}; // struct executor

bool regex_exec(
	const regex_prog & prog, regex_expr & expr, const string_view & text)
{
	executor e;
	auto result = e.regexec(prog, expr, text);
	return result;
}

void regerror(char const * s) { out_printf("re error %s\n", s); }

regex_prog & program::compile(const char * pattern)
{
	static std::unordered_map<std::string, regex_prog_ptr> cache;
	if (cache.count(pattern) == 0)
	{
		cache[pattern] = regex_comp(pattern);
	}
	return *cache[pattern];
}

program::program(const char * pattern) { reset(pattern); }

void program::reset(const char * pattern) { compiled = &compile(pattern); }

program::result_iterator::result_iterator(
	const regex_prog & c, const string_view & s)
	: compiled(&c)
	, rest(s)
{
	advance();
}

void program::result_iterator::advance()
{
	// We start searching for a match at the end of the previous match.
	if (regex_exec(*compiled, expressions, rest))
	{
		// A match means the subexpressions are filled in and the first entry
		// is the full match. Advance `rest` to follow the match.
		rest = string_view(expressions.sub[0].end(), rest.end());
	}
	else
	{
		// No next match, reset to end/nothing.
		rest = string_view{};
		expressions = regex_expr{};
	}
}

}} // namespace b2::regex
