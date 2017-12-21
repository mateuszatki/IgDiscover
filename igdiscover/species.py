"""
Species-specific code, such as lists of motifs and regular expressions.

Some refactoring is needed to make this module actually usable for many
species. Right now, it works for - at least - human, rhesus monkey and mouse.
"""
import re
from sqt.dna import amino_acid_regex
from .utils import nt_to_aa


# Regular expressions for CDR3 detection
#
# The idea comes from D’Angelo et al.: The antibody mining toolbox.
# http://dx.doi.org/10.4161/mabs.27105
# The heavy-chain regex was taken directly from there, but the difference
# is that we express everything in terms of amino acids, not nucleotides.
# This simplifies the expressions and makes them more readable.
#
_CDR3_REGEX = {
	# Heavy chain
	'VH': re.compile("""
		[FY] [FHVWY] C
		(?P<cdr3>
			[ADEGIKMNRSTV] .{3,31}
		)
		W[GAV]
		""", re.VERBOSE),

	# Light chain, kappa
	'VK': re.compile("""
		[FSVY] [CFHNVY] [CDFGLSW]
		(?P<cdr3>
			.{4,15}
		)
		[FLV][GRV]
		""", re.VERBOSE),

	# Light chain, lambda
	'VL': re.compile("""
		# the negative lookahead assertion ensures that the rightmost start is found
		[CDY](?![CDY][CFHSY][CFGW])[CFHSY][CFGW]
		(?P<cdr3>
			.{4,15}
		)
		[FS]G
		""", re.VERBOSE)
}

_CDR3_VH_ALTERNATIVE_REGEX = re.compile("""
		C
		(?P<cdr3> . [RK] .{3,30})
		[WF]G.G
""", re.VERBOSE)


def find_cdr3(sequence, chain):
	"""
	Find the CDR3 in the given sequence, assuming it comes from the given chain ('VH', 'VK', 'VL').
	If the chain is not one of 'VH', 'VK', 'VL', return None.

	Return a tuple (start, stop) if found, None otherwise.
	"""
	try:
		regex = _CDR3_REGEX[chain]
	except KeyError:
		return None
	matches = []
	for offset in 0, 1, 2:
		aa = nt_to_aa(sequence[offset:])
		match = regex.search(aa)
		if not match and chain == 'VH':
			match = _CDR3_VH_ALTERNATIVE_REGEX.search(aa)
		if match:
			start, stop = match.span('cdr3')
			matches.append((start * 3 + offset, stop * 3 + offset))
	return min(matches, default=None)


# Matches the start of the CDR3 within the end of a VH sequence
_CDR3START_VH_REGEX = re.compile("""
	[FY] [FHVWY] C
	(?P<cdr3_start>
		[ADEGIKMNRSTV*] | $
	)
	""", re.VERBOSE)


_CDR3START_VH_ALTERNATIVE_REGEX = re.compile("""
	C
	(?P<cdr3_start> . [RK])
	""", re.VERBOSE)


# Matches after the end of the CDR3 within a J sequence
_CDR3END_JH_REGEX = re.compile('W[GAV]')


def v_cdr3_start(sequence, chain):
	"""
	Find the position of the CDR3 start within the end of a
	V sequence.
	"""
	assert chain == 'VH'
	if 'N' in sequence:
		return None
	aa = nt_to_aa(sequence)
	head, tail = aa[:-12], aa[-12:]
	match = _CDR3START_VH_REGEX.search(tail)
	if not match:
		match = _CDR3START_VH_ALTERNATIVE_REGEX.search(tail)
	if not match:
		return None
	return 3 * (len(head) + match.start('cdr3_start'))


def j_cdr3_end(sequence, chain):
	"""
	Find the position of the CDR3 end within a J sequence

	Return a tuple (frame, cdr3_end) where frame is the frameshift
	(0, 1 or 2).
	"""
	assert chain == 'VH'
	if 'N' in sequence:
		return None
	for frame in 0, 1, 2:
		aa = nt_to_aa(sequence[frame:])
		match = _CDR3END_JH_REGEX.search(aa)
		if match:
			return frame, match.start() * 3 + frame
	return None


# When searching for the CDR3, start this many bases to the left of the end of
# the V match.
CDR3_SEARCH_START = 30
