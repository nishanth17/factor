"""C9a worker settlement and overhead reservations on the durable ledger."""

from ..c10a04.accounting import Ledger


class StudyLedger(Ledger):
    def storage_room(self, stage, reservation):
        used = sum(
            row.get("identity", {}).get("output_byte_reservation", 0)
            for row in self.data["calls"]
            if row["phase"] == stage
        )
        if used + reservation > self.limits[stage]["disk_bytes"]:
            raise RuntimeError("stage output reservation exhausted")

    def call_cell(self, stage, name, cell, function):
        """Reserve startup, warmup and end costs before worker admission."""
        overhead = cell["warm_seconds"]
        used = sum(
            row.get("identity", {}).get("warm_seconds", 0)
            for row in self.data["calls"]
            if row["phase"] == stage
        )
        saved = self.output / (name + ".json")
        if not saved.exists() and used + overhead > self.limits[stage]["warm"]:
            raise RuntimeError("warmup/end reservation exhausted")
        disk = 2**20 + cell.get("observer_byte_limit", 0)
        if not saved.exists():
            self.storage_room(stage, disk)
        work = (
            cell["work_limit"]
            + cell["warm_work_limit"]
            + cell.get("observer_work_limit", 0)
        )

        def checked():
            result = function()
            if (
                not isinstance(result, dict)
                or type(result.get("work")) is not int
            ):
                raise ValueError("worker must propagate exact consumed work")
            if not 0 <= result["work"] <= work:
                raise ValueError("worker work exceeds its reservation")
            if (
                sum(
                    result.get(k, 0)
                    for k in ("measured_work", "warm_work", "observer_work")
                )
                != result["work"]
            ):
                raise ValueError("worker work components do not settle")
            for key, limit in (
                ("measured_work", cell["work_limit"]),
                ("warm_work", cell["warm_work_limit"]),
                ("observer_work", cell.get("observer_work_limit", 0)),
            ):
                if (
                    type(result.get(key)) is not int
                    or not 0 <= result[key] <= limit
                ):
                    raise ValueError("worker component exceeds its allowance")
            return result

        identity = dict(
            cell,
            sample=cell["sample"],
            fixture=cell["fixture"],
            cap=cell["seconds"],
            output_byte_reservation=disk,
        )
        return self.call(
            stage,
            name,
            cell["seconds"] + cell["margin"] + overhead,
            work,
            identity,
            checked,
        )

    def call_operation(self, stage, name, seconds, work, identity, function):
        """Charge analysis, generation and replay against finite allowances."""
        identity = dict(identity, output_byte_reservation=2**20)
        if not (self.output / (name + ".json")).exists():
            self.storage_room(stage, 2**20)

        def checked():
            result = function()
            if (
                type(result.get("work")) is not int
                or not 0 <= result["work"] <= work
            ):
                raise ValueError(
                    "operation omitted or exceeded its work charge"
                )
            return result

        return self.call(stage, name, seconds, work, identity, checked)
