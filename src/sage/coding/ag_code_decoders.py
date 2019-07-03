r"""
Decoders and encoders for AG codes

This module defines decoders and encoders for evaluation AG codes and
differential AG codes.

The algorithms for unique decoding of AG codes implemented in Sage are from
[LBO2014]_, [Lee2016]_.

EXAMPLES::

    sage: F.<a> = GF(9)
    sage: A2.<x,y> = AffineSpace(F, 2)
    sage: C = A2.curve(y^3 + y - x^4)
    sage: pts = [(0,a^2),(0,a^6),(1,2),(1,a),(1,a^3),(2,2),(2,a),
    ....: (2,a^3),(a,1),(a,a^7),(a,a^5), (a^2,2),(a^2,a),(a^2,a^3), (a^7, 1),
    ....: (a^7, a^7), (a^7, a^5), (a^5, 1), (a^5, a^7), (a^5, a^5),
    ....: (a^3,1),(a^3,a^7),(a^3,a^5),(a^6,2),(a^6,a),(a^6,a^3)]
    sage: pts = [C(p) for p in pts]
    sage: pls = [p.places()[0] for p in pts]
    sage: Q = C.places_at_infinity()[0]
    sage: O = C([0,0]).places()[0]
    sage: G = -O + 18*Q
    sage: v = vector([1,a,1,a^7,0,1,a^3,1,a^3,0,a^7,a^7,a^3,1,a^3,0,a^7,0,0,a^3,a^6,a^6,0,a^5,a^7,a^5])
    sage: code = codes.EvaluationAGCode(pls, G)
    sage: dec = code.decoder('uniqueK', Q)
    sage: enc = dec.connected_encoder()
    sage: message = vector([a^6,a^5,2,a^3,a^2,a,1,a^7,a^6,a^5,2,a^3,a^2,a,1])
    sage: dec.decode_to_message(v) == message
    True

AUTHORS:

- Kwankyu Lee (2019-03): initial version

"""
#*****************************************************************************
#       Copyright (C) 2019 Kwankyu Lee <kwankyu@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function

from sage.rings.all import PolynomialRing
from sage.rings.infinity import infinity

from sage.modules.all import vector
from sage.matrix.all import matrix

from .linear_code import LinearCode
from .encoder import Encoder
from .decoder import Decoder, DecodingError


class EvaluationAGCodeEncoder(Encoder):
    """
    Encoder of an evaluation AG code

    INPUT:

    - ``code`` -- an evaluation AG code

    - ``decoder`` -- a decoder of the code

    EXAMPLES::

        sage: F.<a> = GF(4)
        sage: P.<x,y> = AffineSpace(F, 2);
        sage: C = Curve(y^2 + y - x^3)
        sage: F = C.function_field()
        sage: pls = F.places()
        sage: p = C([0,0])
        sage: Q, = p.places()
        sage: D = [pl for pl in pls if pl != Q]
        sage: G = 5*Q
        sage: code = codes.EvaluationAGCode(D, G)
        sage: dec = code.decoder('uniqueK', Q)
        sage: enc = dec.connected_encoder()
        sage: enc
        Encoder for [8, 5] evaluation AG code over GF(4)

    Saving the decoder (and the connected encoder) for later examples and tests::

        sage: save(dec, DOT_SAGE + 'temp/dec')

    The saved object is to be removed after final example.
    """
    def __init__(self, code, decoder):
        """
        Initialize.

        TESTS::

            sage: dec = load(DOT_SAGE + 'temp/dec')
            sage: enc = dec.connected_encoder()
            sage: TestSuite(enc).run(skip='_test_pickling')
        """
        if decoder is None:
            raise TypeError('decoder not yet constructed')

        Encoder.__init__(self, code)

        self._encode = decoder._encode
        self._unencode = decoder._unencode

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: dec = load(DOT_SAGE + 'temp/dec')
            sage: enc = dec.connected_encoder()
            sage: {enc: 1}
            {Encoder for [8, 5] evaluation AG code over GF(4): 1}
        """
        return hash((self.code(), self._encode))

    def __eq__(self, other):
        """
        Test equality.

        EXAMPLES::

            sage: F.<a> = GF(4)
            sage: P.<x,y> = AffineSpace(F, 2);
            sage: C = Curve(y^2 + y - x^3)
            sage: F = C.function_field()
            sage: pls = F.places()
            sage: p = C([0,0])
            sage: Q, = p.places()
            sage: D = [pl for pl in pls if pl != Q]
            sage: G = 5*Q
            sage: code = codes.EvaluationAGCode(D, G)
            sage: dec1 = code.decoder('uniqueK', Q)
            sage: enc1 = dec1.connected_encoder()
            sage: dec2 = code.decoder('uniqueK', Q)
            sage: enc2 = dec2.connected_encoder()
            sage: enc1 == enc2
            True
        """
        if not isinstance(other, EvaluationAGCodeEncoder):
            return False

        return self.code() == other.code() and self._encode == other._encode

    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: dec = load(DOT_SAGE + 'temp/dec')
            sage: enc = dec.connected_encoder()
            sage: enc
            Encoder for [8, 5] evaluation AG code over GF(4)
        """
        return "Encoder for {}" .format(self.code())

    def _latex_(self):
        r"""
        Return the latex representation of ``self``.

        EXAMPLES::

            sage: dec = load(DOT_SAGE + 'temp/dec')
            sage: enc = dec.connected_encoder()
            sage: latex(enc)
            \text{Encoder for }[8, 5]\text{ evaluation AG code over }\Bold{F}_{2^{2}}
        """
        return r"\text{{Encoder for }}{}".format(self.code()._latex_())

    def encode(self, message):
        """
        Return the codeword encoded from the message.

        INPUT:

        - ``message`` -- a vector in the message space

        EXAMPLES::

            sage: dec = load(DOT_SAGE + 'temp/dec')
            sage: enc = dec.connected_encoder()
            sage: msg = enc.message_space().random_element()
            sage: codeword = enc.encode(msg)
            sage: enc.unencode(codeword) == msg
            True
        """
        return self._encode(message)

    def unencode_nocheck(self, codeword):
        """
        Return the message unencoded from ``codeword``.

        INPUT:

        - ``codeword`` -- a vector in the code

        EXAMPLES::

            sage: dec = load(DOT_SAGE + 'temp/dec')
            sage: enc = dec.connected_encoder()
            sage: msg = enc.message_space().random_element()
            sage: codeword = enc.encode(msg)
            sage: enc.unencode(codeword) in enc.message_space()  # indirect doctest
            True

        TESTS::

            sage: import os; os.remove(DOT_SAGE + 'temp/dec.sobj')
        """
        return self._unencode(codeword)


class EvaluationAGCodeUniqueDecoder(Decoder):
    """
    Unique decoder for an evaluation AG code with one extra rational place.

    INPUT:

    - ``code`` -- an evaluation AG code defined by evaluation on places `D`

    - ``Q`` -- a rational place not in `D`

    - ``verbose`` -- boolean; if ``True``, verbose information
      is printed

    EXAMPLES::

        sage: F.<a> = GF(4)
        sage: P.<x,y> = AffineSpace(F, 2);
        sage: C = Curve(y^2 + y - x^3)
        sage: F = C.function_field()
        sage: pls = F.places()
        sage: p = C([0,0])
        sage: Q, = p.places()
        sage: D = [pl for pl in pls if pl != Q]
        sage: G = 5*Q
        sage: code = codes.EvaluationAGCode(D, G)
        sage: dec = code.decoder('uniqueK', Q)
        sage: rv = code.random_received_vector(1)
        sage: dec.decode_to_code(rv) in code
        True

    Saving the decoder for later examples and tests::

        sage: save(dec, DOT_SAGE + 'temp/dec')

    The saved object is to be removed after final example.
    """
    _decoder_type = {'always-succeed'}

    def __init__(self, code, Q, verbose=False):
        """
        Initialize.

        TESTS::

            sage: dec = load(DOT_SAGE + 'temp/dec')
            sage: TestSuite(dec).run()
        """
        if Q.degree() != 1:
            raise ValueError("The place Q is not a rational place")
        if Q in code._pls:
            raise ValueError("THe place Q is one of the places defining the code")

        Decoder.__init__(self, code, code.ambient_space(), connected_encoder_name='evaluation')

        circuit = _Evaluation_AG_Code_Decoder_K(code._pls, code._G, Q, verbose=verbose)

        self._encode = circuit.encode
        self._unencode = circuit.unencode
        self._decoder = circuit.decode_unique_K
        self._info = circuit.info

        self._Q = Q

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: dec = load(DOT_SAGE + 'temp/dec')
            sage: {dec: 1}
            {Unique decoder for [8, 5] evaluation AG code over GF(4): 1}
        """
        return hash((self.code(), self._Q))

    def __eq__(self, other):
        """
        Check whether ``other`` equals ``self``.

        EXAMPLES::

            sage: F.<a> = GF(4)
            sage: P.<x,y> = AffineSpace(F, 2);
            sage: C = Curve(y^2 + y - x^3)
            sage: F = C.function_field()
            sage: pls = F.places()
            sage: p = C([0,0])
            sage: Q, = p.places()
            sage: D = [pl for pl in pls if pl != Q]
            sage: c1 = codes.EvaluationAGCode(D, 5*Q)
            sage: c2 = codes.EvaluationAGCode(D, 5*Q)
            sage: c1 == c2
            True
        """
        if not isinstance(other, type(self)):
            return False
        return self.code() == other.code() and self._Q == other._Q

    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: dec = load(DOT_SAGE + 'temp/dec')
            sage: dec
            Unique decoder for [8, 5] evaluation AG code over GF(4)
        """
        return "Unique decoder for %s" % self.code()

    def _latex_(self):
        r"""
        Return the latex representation of ``self``.

        EXAMPLES::

            sage: dec = load(DOT_SAGE + 'temp/dec')
            sage: latex(dec)
            \text{Unique decoder for }[8, 5]\text{ evaluation AG code over }\Bold{F}_{2^{2}}
        """
        return r"\text{{Unique decoder for }}{}".format(self.code()._latex_())

    def connected_encoder(self, *args, **kwargs):
        r"""
        Returns the connected encoder for this decoder.

        INPUT:

        - ``args``, ``kwargs`` -- all additional arguments are forwarded to the
          constructor of the connected encoder

        EXAMPLES::

            sage: dec = load(DOT_SAGE + 'temp/dec')
            sage: dec.connected_encoder()
            Encoder for [8, 5] evaluation AG code over GF(4)
        """
        return self.code().encoder(self._connected_encoder_name, self, *args, **kwargs)

    def decode_to_code_and_message(self, received_vector, verbose=False):
        r"""
        Return the codeword and the message decoded from ``received_vector``.

        INPUT:

        - ``received_vector`` -- a vector in the ambient space of the code

        - ``verbose`` -- boolean; if ``True``, verbose information on the decoding process
          is printed

        EXAMPLES::

            sage: dec = load(DOT_SAGE + 'temp/dec')
            sage: code = dec.code()
            sage: rv = code.random_received_vector(1)
            sage: cw, msg = dec.decode_to_code_and_message(rv)
            sage: (cw - rv).hamming_weight() == 1
            True
        """
        codeword, message = self._decoder(received_vector, verbose=verbose)
        return codeword, message

    def decode_to_message(self, received_vector, verbose=False):
        r"""
        Return the message decoded from ``received_vector``.

        INPUT:

        - ``received_vector`` -- a vector in the ambient space of the code

        - ``verbose`` -- boolean; if ``True``, verbose information on the decoding process
          is printed

        EXAMPLES::

            sage: dec = load(DOT_SAGE + 'temp/dec')
            sage: enc = dec.connected_encoder()
            sage: code = dec.code()
            sage: rv = code.random_received_vector(1)
            sage: msg = dec.decode_to_message(rv)
            sage: cw = enc.encode(msg)
            sage: (cw - rv).hamming_weight() == 1
            True
        """
        codeword, message = self._decoder(received_vector, verbose=verbose)
        return message

    def decode_to_code(self, received_vector, verbose=False):
        r"""
        Return the codeword decoded from ``received_vector``.

        INPUT:

        - ``received_vector`` -- a vector in the ambient space of the code

        - ``verbose`` -- boolean; if ``True``, verbose information on the decoding process
          is printed

        EXAMPLES::

            sage: dec = load(DOT_SAGE + 'temp/dec')
            sage: code = dec.code()
            sage: rv = code.random_received_vector(1)
            sage: cw = dec.decode_to_code(rv)
            sage: (cw - rv).hamming_weight() == 1
            True

        TESTS::

            sage: import os; os.remove(DOT_SAGE + 'temp/dec.sobj')
        """
        codeword, message = self._decoder(received_vector, verbose=verbose)
        return codeword


class DifferentialAGCodeEncoder(Encoder):
    """
    Encoder of a differential AG code

    INPUT:

    - ``code`` -- a differential AG code

    - ``decoder`` -- a decoder of the code

    EXAMPLES::

        sage: F.<a> = GF(4)
        sage: P.<x,y> = AffineSpace(F, 2);
        sage: C = Curve(y^2 + y - x^3)
        sage: F = C.function_field()
        sage: pls = F.places()
        sage: p = C([0,0])
        sage: Q, = p.places()
        sage: D = [pl for pl in pls if pl != Q]
        sage: G = 5*Q
        sage: code = codes.DifferentialAGCode(D, G)
        sage: dec = code.decoder('uniqueK', Q)
        sage: enc = dec.connected_encoder(); enc
        Residue encoder for [8, 3] differential AG code over GF(4)

    Saving the decoder (and the connected encoder) for later examples and tests::

        sage: save(dec, DOT_SAGE + 'temp/dec')

    The saved object is to be removed after final example.
    """
    def __init__(self, code, decoder=None):
        """
        Initialize.

        TESTS::

            sage: dec = load(DOT_SAGE + 'temp/dec')
            sage: enc = dec.connected_encoder()
            sage: TestSuite(enc).run(skip='_test_pickling')
        """
        if decoder is None:
            raise TypeError('decoder not yet constructed')

        Encoder.__init__(self, code)

        self._encode = decoder._encode
        self._unencode = decoder._unencode

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: dec = load(DOT_SAGE + 'temp/dec')
            sage: enc = dec.connected_encoder()
            sage: {enc: 1}
            {Residue encoder for [8, 3] differential AG code over GF(4): 1}
        """
        return hash((self.code(), self._encode))

    def __eq__(self, other):
        """
        Test equality.

        EXAMPLES::

            sage: F.<a> = GF(4)
            sage: P.<x,y> = AffineSpace(F, 2);
            sage: C = Curve(y^2 + y - x^3)
            sage: F = C.function_field()
            sage: pls = F.places()
            sage: p = C([0,0])
            sage: Q, = p.places()
            sage: D = [pl for pl in pls if pl != Q]
            sage: G = 5*Q
            sage: code = codes.DifferentialAGCode(D, G)
            sage: dec1 = code.decoder('uniqueK', Q)
            sage: enc1 = dec1.connected_encoder()
            sage: dec2 = code.decoder('uniqueK', Q)
            sage: enc2 = dec2.connected_encoder()
            sage: enc1 == enc2
            True
        """
        if not isinstance(other, DifferentialAGCodeEncoder):
            return False

        return self.code() == other.code() and self._encode == other._encode

    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: dec = load(DOT_SAGE + 'temp/dec')
            sage: enc = dec.connected_encoder()
            sage: enc
            Residue encoder for [8, 3] differential AG code over GF(4)
        """
        return "Residue encoder for %s" % self.code()

    def _latex_(self):
        r"""
        Return the latex representation of ``self``.

        EXAMPLES::

            sage: dec = load(DOT_SAGE + 'temp/dec')
            sage: enc = dec.connected_encoder()
            sage: latex(enc)
            \text{Residue encoder for }[8, 3]\text{ differential AG code over }\Bold{F}_{2^{2}}
        """
        return r"\text{{Residue encoder for }}{}".format(self.code()._latex_())

    def encode(self, message):
        """
        Return the codeword encoded from the message.

        INPUT:

        - ``message`` -- a vector in the message space

        EXAMPLES::

            sage: dec = load(DOT_SAGE + 'temp/dec')
            sage: enc = dec.connected_encoder()
            sage: msg = enc.message_space().random_element()
            sage: codeword = enc.encode(msg)
            sage: enc.unencode(codeword) == msg
            True
        """
        return self._encode(message)

    def unencode_nocheck(self, codeword):
        """
        Return the message unencoded from ``codeword``.

        INPUT:

        - ``codeword`` -- a vector in the code

        EXAMPLES::

            sage: dec = load(DOT_SAGE + 'temp/dec')
            sage: enc = dec.connected_encoder()
            sage: msg = enc.message_space().random_element()
            sage: codeword = enc.encode(msg)
            sage: enc.unencode(codeword) in enc.message_space()  # indirect doctest
            True

        TESTS::

            sage: import os; os.remove(DOT_SAGE + 'temp/dec.sobj')
        """
        return self._unencode(codeword)


class DifferentialAGCodeUniqueDecoder(Decoder):
    """
    Unique decoder for a differential AG code with one extra rational place.

    INPUT:

    - ``code`` -- an evaluation AG code defined on places `D`

    - ``Q`` -- a rational place not in `D`

    - ``verbose`` -- boolean (default: ``False``); if ``True``, verbose information
      is printed

    EXAMPLES::

        sage: F.<a> = GF(4)
        sage: P.<x,y> = AffineSpace(F, 2);
        sage: C = Curve(y^2 + y - x^3)
        sage: F = C.function_field()
        sage: pls = F.places()
        sage: p = C([0,0])
        sage: Q, = p.places()
        sage: D = [pl for pl in pls if pl != Q]
        sage: G = 5*Q
        sage: code = codes.DifferentialAGCode(D, G)
        sage: dec = code.decoder('uniqueK', Q)
        sage: rv = code.random_received_vector(2)
        sage: dec.decode_to_code(rv) in code
        True

    Saving the decoder for later examples and tests::

        sage: save(dec, DOT_SAGE + 'temp/dec')

    The saved object is to be removed after final example.
    """
    _decoder_type = {'always-succeed'}

    def __init__(self, code, Q, verbose=False):
        """
        Initialize.

        TESTS::

            sage: dec = load(DOT_SAGE + 'temp/dec')
            sage: TestSuite(dec).run()
        """
        if Q.degree() != 1:
            raise ValueError("The place Q is not a rational place")
        if Q in code._pls:
            raise ValueError("THe place Q is one of the places defining the code")

        Decoder.__init__(self, code, code.ambient_space(), connected_encoder_name='residue')

        circuit = _Differential_AG_Code_Decoder_K(code._pls, code._G, Q, verbose=verbose)

        self._encode = circuit.encode
        self._unencode = circuit.unencode
        self._decoder = circuit.decode_unique_K
        self._info = circuit.info

        self._Q = Q

    def __hash__(self):
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: dec = load(DOT_SAGE + 'temp/dec')
            sage: {dec: 1}
            {Unique decoder for [8, 3] differential AG code over GF(4): 1}
        """
        return hash((self.code(), self._Q))

    def __eq__(self, other):
        """
        Check whether ``other`` equals ``self``.

        EXAMPLES::

            sage: F.<a> = GF(4)
            sage: P.<x,y> = AffineSpace(F, 2);
            sage: C = Curve(y^2 + y - x^3)
            sage: F = C.function_field()
            sage: pls = F.places()
            sage: p = C([0,0])
            sage: Q, = p.places()
            sage: D = [pl for pl in pls if pl != Q]
            sage: c1 = codes.DifferentialAGCode(D, 5*Q)
            sage: c2 = codes.DifferentialAGCode(D, 5*Q)
            sage: c1 == c2
            True
        """
        if not isinstance(other, type(self)):
            return False
        return self.code() == other.code() and self._Q == other._Q

    def _repr_(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES::

            sage: dec = load(DOT_SAGE + 'temp/dec')
            sage: dec
            Unique decoder for [8, 3] differential AG code over GF(4)
        """
        return "Unique decoder for %s" % self.code()

    def _latex_(self):
        r"""
        Return the latex representation of ``self``.

        EXAMPLES::

            sage: dec = load(DOT_SAGE + 'temp/dec')
            sage: latex(dec)
            \text{Unique decoder for }[8, 3]\text{ differential AG code over }\Bold{F}_{2^{2}}
        """
        return r"\text{{Unique decoder for }}{}".format(self.code()._latex_())

    def connected_encoder(self, *args, **kwargs):
        r"""
        Returns the connected encoder for this decoder.

        INPUT:

        - ``args``, ``kwargs`` -- all additional arguments are forwarded to the
          constructor of the connected encoder

        EXAMPLES::

            sage: dec = load(DOT_SAGE + 'temp/dec')
            sage: dec.connected_encoder()
            Residue encoder for [8, 3] differential AG code over GF(4)
        """
        return self.code().encoder(self._connected_encoder_name, self, *args, **kwargs)

    def decode_to_code_and_message(self, received_vector, verbose=False):
        r"""
        Return the codeword and the message decoded from ``received_vector``.

        INPUT:

        - ``received_vector`` -- a vector in the ambient space of the code

        - ``verbose`` -- boolean; if ``True``, verbose information on the decoding process
          is printed

        EXAMPLES::

            sage: dec = load(DOT_SAGE + 'temp/dec')
            sage: code = dec.code()
            sage: rv = code.random_received_vector(2)
            sage: cw, msg = dec.decode_to_code_and_message(rv)
            sage: (cw - rv).hamming_weight() == 2
            True
        """
        codeword, message = self._decoder(received_vector, verbose=verbose)
        return codeword, message

    def decode_to_message(self, received_vector, verbose=False):
        r"""
        Return the message decoded from ``received_vector``.

        INPUT:

        - ``received_vector`` -- a vector in the ambient space of the code

        - ``verbose`` -- boolean; if ``True``, verbose information on the decoding process
          is printed

        EXAMPLES::

            sage: dec = load(DOT_SAGE + 'temp/dec')
            sage: enc = dec.connected_encoder()
            sage: code = dec.code()
            sage: rv = code.random_received_vector(2)
            sage: msg = dec.decode_to_message(rv)
            sage: cw = enc.encode(msg)
            sage: (cw - rv).hamming_weight() == 2
            True
        """
        codeword, message = self._decoder(received_vector, verbose=verbose)
        return message

    def decode_to_code(self, received_vector, verbose=False):
        r"""
        Return the codeword decoded from ``received_vector``.

        INPUT:

        - ``received_vector`` -- a vector in the ambient space of the code

        - ``verbose`` -- boolean; if ``True``, verbose information on the decoding process
          is printed

        EXAMPLES::

            sage: dec = load(DOT_SAGE + 'temp/dec')
            sage: enc = dec.connected_encoder()
            sage: code = dec.code()
            sage: rv = code.random_received_vector(2)
            sage: cw = dec.decode_to_code(rv)
            sage: (cw - rv).hamming_weight() == 2
            True

        TESTS::

            sage: import os; os.remove(DOT_SAGE + 'temp/dec.sobj')
        """
        codeword, message = self._decoder(received_vector, verbose=verbose)
        return codeword


class _Decoder_K(object):
    """
    Common base class for both differential and evaluation AG code decoder K.
    """
    def encode(self, message):
        """
        Encode ``message``.

        TESTS::

            sage: F.<a> = GF(4)
            sage: P.<x,y> = AffineSpace(F, 2);
            sage: C = Curve(y^2 + y - x^3)
            sage: pls = C.places()
            sage: p = C([0,0])
            sage: Q, = p.places()
            sage: D = [pl for pl in pls if pl != Q]
            sage: G = 5*Q
            sage: from sage.coding.ag_code_decoders import _Evaluation_AG_Code_Decoder_K
            sage: circuit = _Evaluation_AG_Code_Decoder_K(D, G, Q)
            sage: rv = vector([0, 0, 0, a, 0, a, a + 1, 0])
            sage: cw, msg = circuit.decode_unique_K(rv)
            sage: circuit.encode(msg) == cw
            True
        """
        code_basis = self.code_basis
        message_index = self.message_index
        return vector(sum([message[i]*code_basis[i] for i in range(len(message_index))]))

    def unencode(self, codeword):
        """
        Unencode ``codeword``.

        TESTS::

            sage: F.<a> = GF(4)
            sage: P.<x,y> = AffineSpace(F, 2);
            sage: C = Curve(y^2 + y - x^3)
            sage: pls = C.places()
            sage: p = C([0,0])
            sage: Q, = p.places()
            sage: D = [pl for pl in pls if pl != Q]
            sage: G = 5*Q
            sage: from sage.coding.ag_code_decoders import _Evaluation_AG_Code_Decoder_K
            sage: circuit = _Evaluation_AG_Code_Decoder_K(D, G, Q)
            sage: rv = vector([0, 0, 0, a, 0, a, a + 1, 0])
            sage: cw, msg = circuit.decode_unique_K(rv)
            sage: circuit.unencode(cw) == msg
            True
        """
        _, message = self.decode_unique_K(codeword)
        return message

    def degree(self, f):
        """
        Return the degree of polynomial ``f``

        TESTS::

            sage: F.<a> = GF(4)
            sage: P.<x,y> = AffineSpace(F, 2);
            sage: C = Curve(y^2 + y - x^3)
            sage: F = C.function_field()
            sage: pls = F.places()
            sage: p = C([0,0])
            sage: Q, = p.places()
            sage: D = [pl for pl in pls if pl != Q]
            sage: G = 5*Q
            sage: from sage.coding.ag_code_decoders import _Evaluation_AG_Code_Decoder_K
            sage: circuit = _Evaluation_AG_Code_Decoder_K(D, G, Q)
            sage: circuit.degree(0)
            -Infinity
        """
        if f.is_zero():
            return -infinity
        else:
            return f.degree()

    def exponents(self, s):
        """
        Return the exponents of the monomial with weighted degree s.

        TESTS::

            sage: F.<a> = GF(4)
            sage: P.<x,y> = AffineSpace(F, 2);
            sage: C = Curve(y^2 + y - x^3)
            sage: F = C.function_field()
            sage: pls = F.places()
            sage: p = C([0,0])
            sage: Q, = p.places()
            sage: D = [pl for pl in pls if pl != Q]
            sage: G = 5*Q
            sage: from sage.coding.ag_code_decoders import _Evaluation_AG_Code_Decoder_K
            sage: circuit = _Evaluation_AG_Code_Decoder_K(D, G, Q)
            sage: circuit.exponents(11)
            (8, 1)
        """
        gamma = self.gamma
        dRbar = self.dRbar  # dWbar for differential AG code
        r = s % gamma
        if s < dRbar[r]:  # if s is a gap, first exponent is set to -1
            return (-1, r)
        else:
            return ((s - dRbar[r]) // gamma, r)

    def substitution(self, vec, w, k, i):
        """
        Substitute z with (z + w*phi_s).

        TESTS::

            sage: F.<a> = GF(4)
            sage: P.<x,y> = AffineSpace(F, 2);
            sage: C = Curve(y^2 + y - x^3)
            sage: F = C.function_field()
            sage: pls = F.places()
            sage: p = C([0,0])
            sage: Q, = p.places()
            sage: D = [pl for pl in pls if pl != Q]
            sage: G = 5*Q
            sage: from sage.coding.ag_code_decoders import _Evaluation_AG_Code_Decoder_K
            sage: circuit = _Evaluation_AG_Code_Decoder_K(D, G, Q)
            sage: W.<x> = F[]
            sage: circuit.substitution(vector([0, a*x^2 + a*x, x + 1, 0]), a, 1, 1)
            (0, a*x^2 + (a*x + a)*x + a*x, x + 1, 0)
        """
        gamma = self.gamma
        mul_mat = self.mul_mat
        x = self.x

        seq = vec.list()
        b = vector(seq[:gamma])
        a = vector(seq[gamma:2*gamma])

        c = (w*x**k)*sum([a[j]*mul_mat[j][i] for j in range(gamma)])
        return vector((c+b).list() + a.list())

    def decode_unique_K(self, received_vector, verbose=False):
        """
        Return the codeword corrected from the received vector.

        INPUT:

        - ``received_vector`` -- a received vector in the ambient space of the
          code

        - ``verbose`` -- boolean; if ``True``, verbose information is printed
          while decoding

        If decoding fails for some reason, ``DecodingError`` is raised. The
        message contained in the exception indicates the type of the decoding
        failure.

        TESTS::

            sage: F.<a> = GF(4)
            sage: P.<x,y> = AffineSpace(F, 2);
            sage: C = Curve(y^2 + y - x^3)
            sage: pls = C.places()
            sage: p = C([0,0])
            sage: Q, = p.places()
            sage: D = [pl for pl in pls if pl != Q]
            sage: G = 5*Q
            sage: from sage.coding.ag_code_decoders import _Evaluation_AG_Code_Decoder_K
            sage: circuit = _Evaluation_AG_Code_Decoder_K(D, G, Q)
            sage: rv = vector(F, [1, a, 1, a + 1, a + 1, a + 1, 1, a + 1])
            sage: circuit.decode_unique_K(rv)
            ((1, a, 1, a + 1, a + 1, 0, 1, a + 1), (1, 0, a + 1, a + 1, a))
        """
        code_length = self.code_length
        designed_distance = self.designed_distance

        gamma = self.gamma
        dR = self.dR
        dRbar = self.dRbar  # dWbar for differential AG code

        hvecs = self.hvecs
        eta_vecs = self.eta_vecs
        mul_mat = self.mul_mat
        coeff_mat = self.coeff_mat

        message_index = self.message_index
        code_basis = self.code_basis

        s0 = self.s0
        tau = self.tau

        W = self.W
        x = self.x

        # auxiliary functions
        degree = self.degree
        exponents = self.exponents
        substitution = self.substitution

        if verbose:
            # auxiliary function for verbose printing
            def vprint_g(g, s):
                if verbose > 1:
                    print(g)
                else:
                    print('[', end='')
                    for i in range(gamma):
                        t = g[i]
                        wd = gamma*degree(t) + dRbar[i]
                        s1 = '{}w{}' if self.is_differential else '{}Y{}'
                        s1 = s1.format(0 if t == 0 else t.lt(), i)
                        if t != 0:
                            s2 = '({})'.format(wd)
                        else:
                            s2 = '(-inf)'
                        print("{:>20} ".format(s1 + s2), end='')
                    for i in range(gamma):
                        t = g[gamma+i]
                        wd = gamma*degree(t) + dR[i] + s
                        s1 = '{}y{}z'.format(0 if t == 0 else t.lt(), i)
                        if t != 0:
                            s2 = '({})'.format(wd)
                        else:
                            s2 = '(-inf)'
                        print('{:>20} '.format(s1 + s2), end='')
                    print(']')

        # construct the initial generators of the interpolation module
        hvec = sum([received_vector[i]*hvecs[i] for i in range(code_length)])

        # weighted degree of hvec
        wd_hvec = max([gamma*degree(hvec[i]) + dRbar[i] for i in range(gamma)])

        message = []
        if wd_hvec <= 0:
            if verbose:
                print("No error.")
            s = 0
            for s in message_index:
                e = exponents(s)
                message.append(hvec[e[1]][e[0]])
        else:
            s = wd_hvec

            mat =[]
            for j in range(gamma):
                row = vector(eta_vecs[j].list() + [W.zero() for i in range(gamma)])
                mat.append(row)
            for j in range(gamma):
                std = [W.zero() for i in range(gamma)]
                std[j] = W.one()
                row = vector(sum([-hvec[i]*mul_mat[j][i] for i in range(gamma)]).list() + std)
                mat.append(row)

            gbmat = matrix(mat)

            nu = []
            for i in range(gamma):
                nu.append(gbmat[i,i].lc())

            found_Q = False
            while s >= s0:
                if verbose:
                    print("s = {}:".format(s))
                    print("Generators:")
                    for j in range(gamma):
                        g = gbmat[j]
                        print("G{} ".format(j), end='')
                        vprint_g(g,s)
                    for j in range(gamma):
                        g = gbmat[gamma+j]
                        print("F{} ".format(j), end='')
                        vprint_g(g,s)

                sk, si = exponents(s)
                delta = 0
                mu = []
                i_k = []
                i_prime = []
                i_value = []
                i_count = []
                voting_value = []
                voting_count = []
                for i in range(gamma):
                    dlt = degree(gbmat[gamma+i,gamma+i])
                    delta += dlt
                    if delta > tau: # more errors than tau; declare failure
                        if verbose:
                            print("decoding failure detected at voting")
                        raise DecodingError("failed at voting")
                    wlt = gamma*dlt + dR[i]
                    if wlt + s + tau < designed_distance:
                        found_Q = True
                        posQ = (s,i)
                        break
                    k,ip = exponents(wlt+s)
                    count = degree(gbmat[ip,ip]) - k
                    i_k.append(k)
                    i_prime.append(ip)
                    i_count.append(count)

                if found_Q:
                    break

                if not s in message_index:
                    for i in range(gamma):
                        k = i_k[i]
                        ip = i_prime[i]

                        if k < 0:
                            value = 0
                        else:
                            value = -gbmat[gamma+i,ip][k]
                        mu.append(1)
                        i_value.append(value)
                    winner = 0
                else: # s is in message_index
                    for i in range(gamma):
                        k = i_k[i]
                        ip = i_prime[i]

                        mui = gbmat[gamma+i,gamma+i].lc()*coeff_mat[i,si]
                        value = -gbmat[gamma+i,ip][k]/mui

                        mu.append(mui)
                        i_value.append(value)

                        cbar = max(i_count[i],0)
                        try:
                            pos = voting_value.index(value)
                            voting_count[pos] += cbar
                        except ValueError:
                            voting_value.append(value)
                            voting_count.append(cbar)

                    # voting
                    c = -1
                    for i in range(len(voting_value)):
                        if c < voting_count[i]:
                            c = voting_count[i]
                            winner = voting_value[i]

                if verbose:
                    print("i_prime:", i_prime)
                    print("i_count:", i_count)
                    print("i_value:", i_value)
                    print("winner = ", winner)

                for i in range(gamma):
                    row_i = gbmat[gamma+i]
                    row_ip = gbmat[i_prime[i]]
                    if winner != 0:
                        row_i = substitution(row_i, winner, sk, si)
                        row_ip = substitution(row_ip, winner, sk, si)
                    if i_value[i] == winner:
                        nrow_ip = row_ip
                        nrow_i = row_i
                    else:
                        nnu = mu[i]*(winner - i_value[i])
                        if i_count[i] > 0:
                            nrow_ip = row_i
                            nrow_i = x**i_count[i]*row_i - nnu/nu[i_prime[i]]*row_ip
                            nu[i_prime[i]] = nnu
                        else:
                            nrow_ip = row_ip
                            nrow_i = row_i - nnu/nu[i_prime[i]]*x**(-i_count[i])*row_ip
                    gbmat[i_prime[i]] = nrow_ip
                    gbmat[gamma+i] = nrow_i

                if s in message_index: # s <= 0 and sk >= 0
                    message.insert(0, winner)

                s -= 1

            if found_Q:
                s,i = posQ
                if verbose:
                    print("found a Q-polynomial at s={}, i={}".format(s,i))
                dlt = gamma*degree(gbmat[gamma+i,gamma+i]) + dR[i]
                while s >= s0:
                    sk, si = exponents(s)
                    if s in message_index: # s le 0 and sk ge 0
                        k, ip = exponents(dlt + s)
                        mui = gbmat[gamma+i,gamma+i].lc()*coeff_mat[i,si]
                        value = -gbmat[gamma+i,ip][k]/mui
                        if not value.is_zero():
                            gbmat[gamma+i] = substitution(gbmat[gamma+i], value, sk, si)
                        message.insert(0,value)
                    s -= 1
                for j in range(gamma):
                    if not gbmat[gamma+i,j].is_zero():
                        if verbose:
                            print("decoding failure detected at division")
                        raise DecodingError("failed at division")

        message = vector(message)
        corrected_vector = vector(sum([message[i]*code_basis[i]
                                       for i in range(len(message_index))]))

        return corrected_vector, message


class _Evaluation_AG_Code_Decoder_K(_Decoder_K):
    """
    Unique decoding algorithm K for evaluation AG codes.

    INPUT:

    - ``pls`` -- list of places of a function field

    - ``G`` -- a divisor of the function field

    - ``Q`` -- a rational place not in ``pls``

    - ``verbose`` -- if ``True``, verbose information is printed.

    TESTS::

        sage: F.<a> = GF(4)
        sage: P.<x,y> = AffineSpace(F, 2);
        sage: C = Curve(y^2 + y - x^3)
        sage: pls = C.places()
        sage: p = C([0,0])
        sage: Q, = p.places()
        sage: D = [pl for pl in pls if pl != Q]
        sage: G = 5*Q
        sage: from sage.coding.ag_code_decoders import _Evaluation_AG_Code_Decoder_K
        sage: circuit = _Evaluation_AG_Code_Decoder_K(D, G, Q)
        sage: rv = vector([a, 0, 0, a, 1, 1, a + 1, 0])
        sage: cw, msg = circuit.decode_unique_K(rv)
        sage: circuit.encode(msg) == cw
        True
    """
    def __init__(self, pls, G, Q, verbose=False):
        """
        Initialize.

        TESTS::

            sage: F.<a> = GF(4)
            sage: P.<x,y> = AffineSpace(F, 2);
            sage: C = Curve(y^2 + y - x^3)
            sage: pls = C.places()
            sage: p = C([0,0])
            sage: Q, = p.places()
            sage: D = [pl for pl in pls if pl != Q]
            sage: G = 5*Q
            sage: from sage.coding.ag_code_decoders import _Evaluation_AG_Code_Decoder_K
            sage: circuit = _Evaluation_AG_Code_Decoder_K(D, G, Q)
            sage: rv = vector([0, 0, 0, a, 0, a, a + 1, 0])
            sage: cw, msg = circuit.decode_unique_K(rv)
            sage: circuit.unencode(cw) == msg
            True
        """
        D = sum(pls)
        F = D.parent().function_field()
        K = F.constant_base_field()
        W = PolynomialRing(K, name='x')  # working polynomial ring
        x = W.gen()

        # length of the code
        code_length = len(pls)

        # compute gamma
        gamma = 1
        while True:
            if Q.divisor(gamma).dimension() > 1:
                break
            gamma += 1

        # compute xR
        for xR in Q.divisor(gamma).basis_function_space():
            if xR.valuation(Q) == -gamma:
                break

        # apery R
        dR = [0 for i in range(gamma)]
        yR = [None for i in range(gamma)]
        s = 0
        n = 0
        while n < gamma:
            g = 0
            for b in Q.divisor(s).basis_function_space():
                if b.valuation(Q) == -s:
                    g = b
                    break
            r = s % gamma
            if g != 0 and not yR[r]:
                dR[r] = s
                yR[r] = g
                n += 1
            s += 1

        # gaps of L
        gaps = set()
        for d in dR:
            gaps.update([d - gamma*(i+1) for i in range(d // gamma)])

        # genus of L
        genus = len(gaps)

        # apery Rbar
        dRbar = [0 for i in range(gamma)]
        yRbar = [None for i in range(gamma)]
        s = -G.degree()
        n = 0
        while n < gamma:
            B = (Q.divisor(s) + G).basis_function_space()
            g = 0
            for b in B:
                if b.valuation(Q) + G.multiplicity(Q) == -s:
                    g = b
                    break
            r = s % gamma
            if g != 0 and not yRbar[r]:
                dRbar[r] = s
                yRbar[r] = g
                n += 1
            s += 1

        if verbose:
            print("gamma = {}".format(gamma))
            print("x = {}".format(xR))
            print("Apery system of R:")
            for i in range(gamma):
                print(" {}: {}, y{} = {}".format(i, dR[i], i, yR[i]))
            print("Apery system of Rbar")
            for i in range(gamma):
                print(" {}: {}, Y{} = {}".format(i, dRbar[i], i, yRbar[i]))

        # ev map for the monomial whose weighted degree is s
        evxR = vector(K, [xR.evaluate(p) for p in pls])
        evyRbar = [vector(K, [yRbar[i].evaluate(p) for p in pls]) for i in range(gamma)]

        self.is_differential = False
        self.code_length = code_length
        self.designed_distance = code_length - G.degree()
        self.gamma = gamma
        self.dR = dR
        self.dRbar = dRbar
        self.W = W
        self.x = x

        degree = self.degree
        exponents = self.exponents

        def monomial(s):
            e = exponents(s)
            return xR**e[0]*yRbar[e[1]]

        def ev_mon(s):
            e = exponents(s)
            return vector([evxR[i]**e[0]*evyRbar[e[1]][i] for i in range(code_length)])

        def next(s):  # integer next to s in Rbar
            while True:
                s += 1
                if exponents(s)[0] >= 0:
                    return s

        # minimum of nongaps of Rbar
        s0 = next(-G.degree() - 1)

        # basis of the code ev(L(G))
        message_index = []
        code_basis = []
        s = s0
        v = ev_mon(s)
        V = v.parent()
        while s <= 0:
            if not V.are_linearly_dependent(code_basis + [v]):
                message_index.append(s)
                code_basis.append(v)
            s = next(s)
            v = ev_mon(s)

        ## underlying linear code, not used anywhere
        #self.code = LinearCode(matrix(code_basis))

        # compute a basis of J and h-functions via FGLM algorithm
        def get_eta_basis():
            basis = [None for i in range(gamma)]
            s = s0
            mat = matrix(ev_mon(s))
            delta = [exponents(s)]
            num = 0
            while num < gamma:
                s = next(s)
                e = exponents(s)
                if basis[e[1]] is None:
                    v = ev_mon(s)
                    try:
                        sol = mat.solve_left(v)
                        gen = [W.zero() for i in range(gamma)]
                        for i in range(len(delta)):
                            gen[delta[i][1]] += -sol[i]*x**delta[i][0]
                        gen[e[1]] += x**e[0]
                        basis[e[1]] = vector(gen)
                        num += 1
                    except ValueError:
                        mat = matrix(list(mat) + [v])
                        delta.append(e)

            vecs = [None for i in range(code_length)]
            matinv = mat.inverse()
            for i in range(code_length):
                h = [W.zero() for k in range(gamma)]
                for j in range(code_length):
                    h[delta[j][1]] += matinv[i,j]*x**delta[j][0]
                vecs[i] = vector(h)

            return basis, vecs

        eta_vecs, hvecs = get_eta_basis()

        if verbose:
            print("message indices:")
            print(message_index)
            print("eta basis:")
            print(eta_vecs)
            print("Lagrange polynomials:")
            for i in range(code_length):
                print("h{} = {}".format(i, hvecs[i]))

        # Lee-O'Sullivan bound
        def nu(s):
            m = 0
            for i in range(gamma):
                e = exponents(s + dR[i])
                m += max(0, degree(eta_vecs[e[1]][e[1]]) - e[0])
            return m

        nus = [nu(s) for s in message_index]
        dLO = min(nus)
        tau = (dLO - 1) // 2

        if verbose:
            print("nu's =", nus)
            print("d_LO =", dLO)
            print("tau =", tau)

        # the vector form corresponding to f in Rbar
        def vec_form(f):
            r = f
            l = [W.zero() for i in range(gamma)]
            while r != 0:
                s = -r.valuation(Q) - G.multiplicity(Q)
                e = exponents(s)
                mon = xR**e[0]*yRbar[e[1]]
                c = (r/mon).evaluate(Q)
                l[e[1]] += c*x**e[0]
                r -= c*mon
            return vector(l)

        # the matrix of the leading coefficient of y_i*ybar_j and the product
        def get_mul_mat():
            cm = matrix.zero(K, gamma, gamma)
            vm = [[None for j in range(gamma)] for i in range(gamma)]
            for i in range(gamma):
                for j in range(gamma):
                    f = yR[i]*yRbar[j]
                    v = vec_form(f)
                    e = exponents(dR[i]+dRbar[j])
                    cm[i,j] = v[e[1]][e[0]]
                    vm[i][j] = v
            return vm, cm

        mul_mat, coeff_mat = get_mul_mat()

        if verbose:
            print("multiplication table:")
            for i in range(gamma):
                for j in range(gamma):
                    print("y{} * Y{}:".format(i, j), mul_mat[i][j])
            print("coefficient array:")
            print(coeff_mat)

        self.code_basis = code_basis
        self.message_index = message_index
        self.hvecs = hvecs
        self.eta_vecs = eta_vecs
        self.mul_mat = mul_mat
        self.coeff_mat = coeff_mat
        self.s0 = s0
        self.tau = tau

        info = {}
        info['gaps'] = gaps
        info['genus'] = genus
        info['nus'] = nus
        info['designed_distance'] = dLO
        info['decoding_radius'] = tau

        self.info = info


class _Differential_AG_Code_Decoder_K(_Decoder_K):
    """
    Unique decoding algorithm K for differential AG codes.

    INPUT:

    - ``pls`` -- list of places of a function field

    - ``G`` -- a divisor of the function field

    - ``Q`` -- a rational place not in ``pls``

    - ``verbose`` -- if ``True``, verbose information is printed.

    TESTS::

        sage: F.<a> = GF(4)
        sage: P.<x,y> = AffineSpace(F, 2);
        sage: C = Curve(y^2 + y - x^3)
        sage: pls = C.places()
        sage: p = C([0,0])
        sage: Q, = p.places()
        sage: D = [pl for pl in pls if pl != Q]
        sage: G = 5*Q
        sage: from sage.coding.ag_code_decoders import _Differential_AG_Code_Decoder_K
        sage: circuit = _Differential_AG_Code_Decoder_K(D, G, Q)
        sage: rv = vector([1, a, 1, a, 1, a, a, a + 1])
        sage: cw, msg = circuit.decode_unique_K(rv)
        sage: circuit.encode(msg) == cw
        True
    """
    def __init__(self, pls, G, Q, verbose=False):
        """
        Initialize

        TESTS::

            sage: F.<a> = GF(4)
            sage: P.<x,y> = AffineSpace(F, 2);
            sage: C = Curve(y^2 + y - x^3)
            sage: pls = C.places()
            sage: p = C([0,0])
            sage: Q, = p.places()
            sage: D = [pl for pl in pls if pl != Q]
            sage: G = 5*Q
            sage: from sage.coding.ag_code_decoders import _Differential_AG_Code_Decoder_K
            sage: circuit = _Differential_AG_Code_Decoder_K(D, G, Q)
            sage: rv = vector([a + 1, a + 1, 0, 1, a, 1, 1, 0])
            sage: cw, msg = circuit.decode_unique_K(rv)
            sage: circuit.unencode(cw) == msg
            True
        """
        D = sum(pls)
        F = D.parent().function_field()
        K = F.constant_base_field()
        W = PolynomialRing(K, name='x')  # working polynomial ring
        x = W.gen()

        # length of the code
        code_length = len(pls)

        # compute gamma
        gamma = 1
        while True:
            if Q.divisor(gamma).dimension() > 1:
                break
            gamma += 1

        # compute xR
        for xR in Q.divisor(gamma).basis_function_space():
            if xR.valuation(Q) == -gamma:
                break

        # apery R
        dR = [0 for i in range(gamma)]
        yR = [None for i in range(gamma)]
        s = 0
        n = 0
        while n < gamma:
            g = 0
            for b in Q.divisor(s).basis_function_space():
                if b.valuation(Q) == -s:
                    g = b
                    break
            r = s % gamma
            if g != 0 and not yR[r]:
                dR[r] = s
                yR[r] = g
                n += 1
            s += 1

        # gaps of L
        gaps = set()
        for d in dR:
            gaps.update([d - gamma*(i+1) for i in range(d // gamma)])

        # genus of L
        genus = len(gaps)

        # apery Wbar
        dWbar = [0 for i in range(gamma)]
        wWbar = [None for i in range(gamma)]
        s = -code_length + G.degree() - 2*genus + 2
        n = 0
        while n < gamma:
            B = (-D + G - Q.divisor(s)).basis_differential_space()
            g = 0
            for b in B:
                if b.valuation(Q) == G.multiplicity(Q) - s:
                    g = b
                    break
            r = s % gamma
            if g != 0 and not wWbar[r]:
                dWbar[r] = s
                wWbar[r] = g
                n += 1
            s += 1

        if verbose:
            print("gamma = {}".format(gamma))
            print("x = {}".format(xR))
            print("Apery system of R")
            for i in range(gamma):
                print(" {}: {}, y{} = {}".format(i, dR[i], i, yR[i]))
            print("Apery system of Wbar")
            for i in range(gamma):
                print(" {}: {}, w{} = {}".format(i, dWbar[i], i, wWbar[i]))

        # res map for the monomial whose weighted degree is s
        evxR = vector(K, [xR.evaluate(p) for p in pls])
        reswWbar = [vector(K, [wWbar[i].residue(p) for p in pls]) for i in range(gamma)]

        self.is_differential = True
        self.code_length = code_length
        self.designed_distance = G.degree() - 2*genus + 2
        self.gamma = gamma
        self.dR = dR
        self.dRbar = dWbar
        self.W = W
        self.x = x

        degree = self.degree
        exponents = self.exponents

        def monomial(s):
            e = exponents(s)
            return xR**e[0]*wWbar[e[1]]

        def res_mon(s):
            e = exponents(s)
            return vector([evxR[i]**e[0]*reswWbar[e[1]][i] for i in range(code_length)])

        def next(s):  # integer next to s in Wbar
            while True:
                s += 1
                if exponents(s)[0] >= 0:
                    return s

        # minimum of nongaps of Wbar
        s0 = next(-code_length + G.degree() - 2*genus + 1)

        # basis of the code res(Omega(G))
        message_index = []
        code_basis = []
        s = s0
        v = res_mon(s)
        V = v.parent()
        while s <= 0:
            if not V.are_linearly_dependent(code_basis + [v]):
                message_index.append(s)
                code_basis.append(v)
            s = next(s)
            v = res_mon(s)

        ## underlying linear code, not used anywhere
        #self.code = LinearCode(matrix(code_basis))

        # compute a basis of J and h-functions via FGLM algorithm
        def get_eta_basis():
            basis = [None for i in range(gamma)]
            s = s0
            mat = matrix(res_mon(s))
            delta = [exponents(s)]
            num = 0
            while num < gamma:
                s = next(s)
                e = exponents(s)
                if basis[e[1]] is None:
                    v = res_mon(s)
                    try:
                        sol = mat.solve_left(v)
                        gen = [W.zero() for i in range(gamma)]
                        for i in range(len(delta)):
                            gen[delta[i][1]] += -sol[i]*x**delta[i][0]
                        gen[e[1]] += x**e[0]
                        basis[e[1]] = vector(gen)
                        num += 1
                    except ValueError:
                        mat = matrix(list(mat) + [v])
                        delta.append(e)

            vecs = [None for i in range(code_length)]
            matinv = mat.inverse()
            for i in range(code_length):
                h = [W.zero() for k in range(gamma)]
                for j in range(code_length):
                    h[delta[j][1]] += matinv[i,j]*x**delta[j][0]
                vecs[i] = vector(h)

            return basis, vecs

        eta_vecs, hvecs = get_eta_basis()

        if verbose:
            print("message indices:")
            print(message_index)
            print("eta basis:")
            print(eta_vecs)
            print("Lagrange polynomials:")
            for i in range(code_length):
                print("h{} = {}".format(i, hvecs[i]))

        # Lee-O'Sullivan bound
        def nu(s):
            m = 0
            for i in range(gamma):
                e = exponents(s + dR[i])
                m += max(0, degree(eta_vecs[e[1]][e[1]]) - e[0])
            return m

        nus = [nu(s) for s in message_index]
        dLO = min(nus)
        tau = (dLO - 1) // 2

        if verbose:
            print("nu's =", [nu(s) for s in message_index])
            print("d_LO =", dLO)
            print("tau =", tau)

        # the vector form corresponding to f in Wbar
        def vec_form(f):
            r = f
            l = [W.zero() for i in range(gamma)]
            while r != 0:
                s = -r.valuation(Q) + G.valuation(Q)
                e = exponents(s)
                mon = xR**e[0]*wWbar[e[1]]
                c = (r/mon).evaluate(Q)
                l[e[1]] += c*x**e[0]
                r -= c*mon
            return vector(l)

        # the matrix of the leading coefficient of y_i*w_j and the product
        def get_mul_mat():
            cm = matrix.zero(K, gamma, gamma)
            vm = [[None for j in range(gamma)] for i in range(gamma)]
            for i in range(gamma):
                for j in range(gamma):
                    f = yR[i]*wWbar[j]
                    v = vec_form(f)
                    e = exponents(dR[i]+dWbar[j])
                    cm[i,j] = v[e[1]][e[0]]
                    vm[i][j] = v
            return vm, cm

        mul_mat, coeff_mat = get_mul_mat()

        if verbose:
            print("multiplication table:")
            for i in range(gamma):
                for j in range(gamma):
                    print("y{} * w{}:".format(i, j), mul_mat[i][j])
            print("coefficient array:")
            print(coeff_mat)

        self.code_basis = code_basis
        self.message_index = message_index
        self.hvecs = hvecs
        self.eta_vecs = eta_vecs
        self.mul_mat = mul_mat
        self.coeff_mat = coeff_mat
        self.s0 = s0
        self.tau = tau

        info = {}
        info['gaps'] = gaps
        info['genus'] = genus
        info['nus'] = nus
        info['designed_distance'] = dLO
        info['decoding_radius'] = tau

        self.info = info
