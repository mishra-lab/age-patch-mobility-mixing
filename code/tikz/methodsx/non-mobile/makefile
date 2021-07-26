latex:
	echo "\halfhometrue" > case.tex && pdflatex main.tex && mv main.pdf pools-half-home.pdf && \
	echo "\halfstaytrue" > case.tex && pdflatex main.tex && mv main.pdf pools-half-stay.pdf && \
	echo "\diffhometrue" > case.tex && pdflatex main.tex && mv main.pdf pools-diff-home.pdf && \
	echo "\diffstaytrue" > case.tex && pdflatex main.tex && mv main.pdf pools-diff-stay.pdf && \
	echo "\halfhometrue"

proportions:
	python3.8 props.py
