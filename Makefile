.PHONY: build bash clean

build:
	docker build -t amplipore .

bash:
	docker run -v ./app:/usr/src/app/ -v ./data:/usr/src/data/ -p 8000:8000 -it --rm amplipore /bin/bash

clean:
	docker system prune --volumes -f