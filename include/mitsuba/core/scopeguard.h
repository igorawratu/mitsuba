#pragma once

#include <utility>

struct ScopeGuardBase {
	ScopeGuardBase() : active_(true){}

	ScopeGuardBase(ScopeGuardBase&& other) : active_(other.active_) {
		other.dismiss();
	}

	void dismiss() noexcept	{
		active_ = false;
	}

protected:
	~ScopeGuardBase() = default;
	bool active_;
};

template<class Func>
struct ScopeGuard : public ScopeGuardBase {
	ScopeGuard() = delete;
	ScopeGuard(const ScopeGuard&) = delete;

	ScopeGuard(Func f) noexcept : 
		ScopeGuardBase(),
		f_(std::move(f)){}

	ScopeGuard(ScopeGuard&& other) noexcept :
		ScopeGuardBase(std::move(other)),
		f_(std::move(other.f_)){}

	~ScopeGuard() noexcept {
		if (active_) {
			try {
				f_();
			}
			catch (...){
			}
		}
	}

	ScopeGuard& operator=(const ScopeGuard&) = delete;

private:
	Func f_;
};

template<class Func> ScopeGuard<Func> makeScopeGuard(Func f) {
	return ScopeGuard<Func>(std::move(f));
}