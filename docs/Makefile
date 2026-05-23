EXEC=rst2html
FILES=emfm math optimize qed_mcpm quantum rochester vdm

all:
	for file in $(FILES); do \
		$(EXEC) $${file}.rst $${file}.html; \
	done

clean:
	for file in $(FILES); do \
		rm -f $${file}.html; \
	done
