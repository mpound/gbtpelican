#!/bin/tcsh -f
foreach session (1 3 4 5 7)
    #echo ./tsys.py --no-tsys --beams --freeze --no-vc3 AGBT25B_386_0${session}
    #./tsys.py --no-tsys --beams --freeze --no-vc3 AGBT25B_386_0${session}
    echo ./tsys.py --parse-only AGBT25B_386_0${session} > tsyslogs/beam_stats_0${session}.md
    ./tsys.py --parse-only AGBT25B_386_0${session} > tsyslogs/beam_stats_0${session}.md
end
