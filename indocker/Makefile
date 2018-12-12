TASK1_PATH=/Users/sokks/Documents/sem7/quantum/prac_scalapack/task1


task1: task1_create_volume task1_run_c

task1_run_c:
	docker container stop task1 || true && docker container rm task1 || true
	docker image rm task1 || true && docker build -t task1 . || true
	docker run -it --name task1 -v task1_vol:/root/task1_vol task1

task1_create_volume:
	docker volume create --driver local \
      --opt type=none \
      --opt device=$(TASK1_PATH) \
      --opt o=bind \
      task1_vol

task1_attach:
	docker attach task1
