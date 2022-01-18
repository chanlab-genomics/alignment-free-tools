import multiprocessing
import time


def worker(num, my_queue):
    name = multiprocessing.current_process().name
    print 'adding', num

    my_queue.put(num)


def my_service(num, my_queue):
    name = multiprocessing.current_process().name
    print 'adding', num
    my_queue.put(num)


if __name__ == '__main__':

    my_queue = multiprocessing.Queue()

    service = multiprocessing.Process(
        name='my_service', target=my_service, args=(1, my_queue))
    worker_1 = multiprocessing.Process(
        name='worker 1', target=worker, args=(2, my_queue))
    worker_2 = multiprocessing.Process(
        target=worker, args=(3, my_queue))  # use default name

    worker_1.start()
    worker_2.start()
    service.start()

    worker_1.join()
    worker_2.join()
    service.join()

    queue_list = [ my_queue.get_nowait() for _ in range(my_queue.qsize())]

    print(queue_list)
    # print('GOT: ', my_queue.get_nowait())
    # print('GOT: ', my_queue.get_nowait())
    # print('GOT: ', my_queue.get_nowait())
    # print(sum(list(my_queue)))
