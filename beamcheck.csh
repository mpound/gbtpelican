#!/bin/tcsh -f
foreach session (1 3 4 5 7)
    echo ./dotsys --no-tsys --beams --freeze --no-vc3 AGBT25B_376_0${session}
end
