"""
Testers for linear codes

EXAMPLES::

    sage: from sage.coding.tester import test_decoder
    sage: F.<a> = GF(4)
    sage: P.<x,y> = AffineSpace(F, 2);
    sage: C = Curve(y^2 + y - x^3)
    sage: F = C.function_field()
    sage: pls = F.places()
    sage: Q, = C([0,0]).places()
    sage: D = [pl for pl in pls if pl != Q]
    sage: code = codes.EvaluationAGCode(D, 3*Q)
    sage: dec = code.decoder('uniqueK', Q)
    sage: dec
    Unique decoder for [8, 3] evaluation AG code over GF(4)

Testing the decoder within the decoding radius::

    sage: test_decoder(dec, 2, 10)
    Success profile: correct 10, incorrect None
    Failure profile: None
    Average decoding time per frame: ...

Testing again out of the decoding radius::

    sage: test_decoder(dec, 3, 10)  # random
    Success profile: correct 2, incorrect {2: 4, 3: 1}
    Failure profile: {'failure at division': 3}
    Average decoding time per frame: 0.0938488006592

Testing with verbose information::

    sage: test_decoder(dec, 4, 10, verbose=True)  # random
    Testing 10 frames with 4 errors:
    #0 incorrect
    #1 incorrect
    #2 incorrect
    #3 failed - failure at division
    #4 incorrect
    #5 incorrect
    #6 incorrect
    #7 failed - failure at division
    #8 incorrect
    #9 failed - failure at division
    Success profile: correct 0, incorrect {2: 4, 3: 2, 4: 1}
    Failure profile: {'failure at division': 3}
    Average decoding time per frame: 0.0874499082565

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
import time
import six
import sys

from collections import Counter

from .decoder import DecodingError

if six.PY2:
    timer = time.time
else:
    timer = time.process_time

def test_decoder(decoder, weight, number, verbose=False):
    """
    Print the result of decoding random received vectors with errors ``weight``
    for ``number`` frames.

    INPUT:

    - ``decoder`` -- decoder of a linear code

    - ``weight`` -- positive integer; weight of random error vector

    - ``number`` -- number of frames to test

    - ``verbose`` -- boolean; if ``True``, verbose information is printed

    EXAMPLES::

        sage: from sage.coding.tester import test_decoder
        sage: F.<a> = GF(4)
        sage: P.<x,y> = AffineSpace(F, 2);
        sage: C = Curve(y^2 + y - x^3)
        sage: F = C.function_field()
        sage: pls = F.places()
        sage: Q, = C([0,0]).places()
        sage: D = [pl for pl in pls if pl != Q]
        sage: code = codes.DifferentialAGCode(D, 3*Q)
        sage: dec = code.decoder('uniqueK', Q)
        sage: dec
        Unique decoder for [8, 5] differential AG code over GF(4)

    ::

        sage: test_decoder(dec, 0, 10)
        Success profile: correct 10, incorrect None
        Failure profile: None
        Average decoding time per frame: ...

    ::

        sage: test_decoder(dec, 1, 20)
        Success profile: correct 20, incorrect None
        Failure profile: None
        Average decoding time per frame: ...

    ::

        sage: test_decoder(dec, 2, 30)  # random
        Success profile: correct 3, incorrect {1: 5, 2: 12}
        Failure profile: {'failed at voting': 10}
        Average decoding time per frame: 0.0731041034063

    The last profiles report that 3 frames were decoded correctly and 3 frames
    were decoded correctly but to codewords of distance `5` or `12` from the
    received vectors respectively, and that the decoder failed to decode 9
    frames declaring ``failed at voting``.
    """
    code = decoder.code()

    total_elapsed_time = 0
    correct = []
    incorrect = []
    failure = []

    if verbose:
        print("Testing {}\nDecoding {} frames with {} errors:".format(decoder, number, weight))

    for frame in range(number):
        if verbose:
            sys.stdout.write('#{} '.format(frame))
        received_vector, sent_codeword = code.random_received_vector_and_codeword(weight)
        try:
            start = timer()
            codeword, message = decoder.decode_to_code_and_message(received_vector)
            end = timer()
            elapsed_time = end - start
            check = (sent_codeword == codeword)
            if check:
                correct.append(True)
                if verbose:
                    print('correct')
            else:
                incorrect.append((received_vector - codeword).hamming_weight())
                if verbose:
                    print('incorrect')
        except DecodingError as e:
            end = timer()
            elapsed_time = end - start
            failure_message = e.args[0]
            failure.append(failure_message)
            if verbose:
                print('failed - {}'.format(failure_message))

        total_elapsed_time += elapsed_time  # elapsed time for decoding

    correct = Counter(correct)
    incorrect = Counter(incorrect)
    failure = Counter(failure)

    correct_tally = correct[True]
    incorrect_tally = dict(incorrect) if len(incorrect) > 0 else None
    failure_tally = dict(failure) if len(failure) > 0 else None

    print("Success profile: correct {}, incorrect {}".format(correct_tally, incorrect_tally))
    print("Failure profile: {}".format(failure_tally))
    print("Average decoding time per frame: {}".format(total_elapsed_time / number))

