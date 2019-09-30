#ifndef BLOCKINGQUEUE_H_
#define BLOCKINGQUEUE_H_

#include <queue>
#include <mutex>
#include <condition_variable>
#include <tuple>

template <class T, class Container = std::deque<T>>
class BlockingQueue{
public:
    BlockQueue(std::uint32_t max_size = 0) : max_size_(max_size), open_(true){}

    ~BlockQueue(){}
    
    void close(){
        {
            std::lock_guard<std::mutex> lock(mutex_);
            open_ = false;
        }
        pull_cond_.notify_all();
    }

    void push(T&& item){
        {
            std::lock_guard<std::mutex> lock(mutex_);
            
            if(max_size_ > 0){
                push_cond_.wait(lock, [this](){return this->queue.size() < max_size_;});
            }

            queue_.push(item);
        }

        pull_cond_.notify_one();
    }

    bool pull(T& item){
        {
            std::lock_guard<std::mutex> lock(mutex_);

            pull_cond_.wait(lock, [this](){return !this->queue_.empty() || !this->open_;});
            if(queue_.empty()){
                return false;
            }
            
            item = std::move(queue_.front());
            queue_.pop_front();
        }

        push_cond_.notify_one();
        return true;
    }

private:
    std::queue<T, Container> queue_;
    std::mutex mutex_;
    std::condition_variable push_cond_;
    std::condition_variable pull_cond_;
    std::uint32_t max_size_;
    bool open_;
}

#endif