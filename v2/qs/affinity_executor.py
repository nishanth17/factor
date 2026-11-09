"""One bounded process per logical slot for current-QS family placement."""


class FamilyAffinityExecutor:
    """Pin packets explicitly; round-robin broadcast calls visit every slot.

    Installation and final CPU snapshots already submit exactly N barrier
    calls. Round-robin placement gives each of those calls a distinct worker.
    Packet placement uses the original family index, leaving packet identities,
    ordered admission, leases and checkpoint coverage unchanged.
    """

    def __init__(
        self, workers, *, mp_context, initializer, initargs, executor_type=None
    ):
        if executor_type is None:
            from .parallel import ObservedProcessPoolExecutor

            executor_type = ObservedProcessPoolExecutor
        self.executors = [
            executor_type(
                max_workers=1,
                mp_context=mp_context,
                initializer=initializer,
                initargs=initargs,
            )
            for _ in range(workers)
        ]
        self.next_slot = 0

    @property
    def resource_owner(self):
        return self.executors[0].resource_owner

    @resource_owner.setter
    def resource_owner(self, budget):
        for executor in self.executors:
            executor.resource_owner = budget

    @property
    def _processes(self):
        # Each manager owns a separate table and can clear it after a crash.
        # CollectionPool's independent creation registry retains lost handles.
        return {
            pid: process
            for executor in self.executors
            for pid, process in (executor._processes or {}).copy().items()
        }

    @property
    def _executor_manager_threads(self):
        return tuple(
            executor._executor_manager_thread
            for executor in self.executors
            if executor._executor_manager_thread is not None
        )

    def submit(self, function, *args, **kwargs):
        slot = self.next_slot
        self.next_slot = (slot + 1) % len(self.executors)
        return self.submit_to(slot, function, *args, **kwargs)

    def submit_to(self, slot, function, *args, **kwargs):
        return self.executors[slot].submit(function, *args, **kwargs)

    def shutdown(self, wait=True, *, cancel_futures=False):
        for executor in self.executors:
            executor.shutdown(wait=wait, cancel_futures=cancel_futures)
