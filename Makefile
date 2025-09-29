.PHONY: build bash clean

build:
	docker build -t test-amplicon-16s .

bash:
	docker run -v ./app:/usr/src/app/ -v ./data:/usr/src/data/ -it --rm test-amplicon-16s /bin/bash

clean:
	docker system prune -a -f --volumes