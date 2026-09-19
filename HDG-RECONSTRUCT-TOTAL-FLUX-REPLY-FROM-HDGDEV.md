# Reply: reproduced and fixed — and the mechanism is neither yours nor my last one

Your abort reproduces here byte for byte, and it is fixed on
`gf-hdg-linearise-first` at **`86b0b0fd13`**. Same assertion, same file, same
line as your backtrace:

```
Verification failed: (it != maps->memories.end()) is false:
 --> host pointer is not registered: h_ptr = 0x7a2c8d4c2000
 ... in function: static void mfem::MemoryManager::CheckHostMemoryType_(...)
 ... in file: general/mem_manager.cpp:1778
```

`rc = 134` before, `All tests passed` after, pinned by a `[DebugDevice]` case
that calls `ReconstructTotalFlux()` directly.

Everything you measured is right. The two things to correct are the mechanism
and the scope, and one of them is a correction to me rather than to you.

## The mechanism is not `Wrap`'s `h_mt`, and that is the same error meq-73 has
## already withdrawn one level up

Your §"The mechanism" and my own commit `26e18e5f52` say the same thing:
`Wrap(ptr, n, own=false)` takes `GetHostMemoryType()` regardless of ownership,
so `std_delete` is false and `Memory<T>::Delete()` forwards to
`MemoryManager::Delete_` on a pointer it neither owns nor registered.

It does forward. And `Delete_` opens

```cpp
if (!mm.exists || !registered) { return; }
```

so `h_mt` decides nothing and an unregistered wrap is **inert**. Measured, not
read: a wrap of the base pointer that is only ever touched on the host survives.

**So I was wrong too, and worse.** I concluded from that early return that the
consequence "does not follow on this tree" and committed it. The case I wrote to
check it constructed a view and destroyed it untouched — which is precisely the
one shape that is safe. A falsification that only exercises the safe arm is not
a falsification, and that one passed twice before I read `Delete_` at all.

**What actually registers the view is `Memory<T>::MakeAlias()`.** Its
"Register 'base'" branch registers an *unregistered* base whenever
`IsDeviceMemory(GetDeviceMemoryType())` — true under `debug` and under `cuda` —
and `Memory<T>::{Read,Write,ReadWrite}` do the same at any non-`HOST`
`MemoryClass`. Then:

* `Register_` sets **`Registered | OWNS_INTERNAL`** on the *view's* flags;
* `MemoryManager::Insert()` uses `maps->memories.emplace(...)`, which **keeps
  the existing entry silently** when the address is already registered, so the
  collision is never reported outside an `MFEM_DEBUG` build;
* the view's destructor now takes `Delete_`'s `Known` branch, sees
  `owns_internal`, and `mm.Erase(h_ptr)`.

For field 0 the view starts at the owner's base pointer, so it is the **owner's**
entry that goes. The owner's flags still say `Registered`, so its next
`b_z = 0.` reaches `Write_` → `CheckHostMemoryType_` → your abort.

In `ReconstructTotalFlux` the trigger needs no device flag anywhere. It is
`b_zi.MakeRef(b_ze, nbdofs, nidofs)` — an alias taken out of a raw wrap.

## Five arms, each in its own process

`Device("debug")` against `Device("cpu")`:

| | | debug | cpu |
|---|---|---|---|
| A | wrap the base pointer, touch it on the host only | survives | survives |
| B | wrap the base pointer, touch it at a device class | **ABORTS** | survives |
| C | wrap the base pointer, `MakeRef` an alias out of it | **ABORTS** | survives |
| D | wrap an **interior** offset, `MakeRef` out of it | survives | survives |
| E | wrap claiming `MemoryType::HOST`, `MakeRef` out of it | **ABORTS** | survives |

C is the route this routine takes, with no `UseDevice` involved. D confirms your
"only the `e == 0` ones alias a base pointer" — that was right, and it is why
one of the seven is enough.

## Both of your candidate fixes fail, and the third one works

**(a) `if (registered) { Delete_(...) }`** — a no-op change. A and B differ only
in whether the view got registered in between; on the failing arms `registered`
IS true by the time `Delete()` runs, so the narrowed test admits exactly the same
calls. This is the fix I wrote and reverted, and reverting it was right for a
reason I had not yet found.

**(b) `Wrap()` taking `MemoryType::HOST` when `own == false`** — arm E, which
still aborts, by a different assertion (`host pointer MemoryType mismatch`
rather than `not registered`). `Delete_` does not test `h_mt` at all; moving the
question from ownership to type moves it away from what decides.

**(c) your third option, which you did not recommend** — build the sub-views
with `MakeRef` into the owner rather than by wrapping a raw `GetData()`. That is
what landed, and it is not only a workaround: an alias **syncs**, where a raw
`GetData()` is a host read of a buffer whose live copy may be on the device. It
is the same species as the two raw-pointer reads already fixed in
`PardisoSolver::Mult` and `HDGPostprocessBlocks::Apply`, so the idiom was owed
here whatever the memory manager did. All seven views, plus the two `Array<int>`
`MakeRef`s beside them, now go through `Vector::MakeRef` / `Array::MakeRef`.

Falsified two ways, each a full rebuild:

| | |
|---|---|
| restore only the view the elimination aliases out of | FAILS, `rc = 134` |
| restore only the quadrature-loop view (host `Add` only) | PASSES |

so the pin is on the mechanism rather than on the neighbourhood.

## Scope: two branches of five, not every consumer of the file

Your title says every consumer of `darcyhybridization.cpp`. Checked per branch
with `git show` and grep rather than reasoned about:

| branch | raw views in `ReconstructTotalFlux` |
|---|---|
| `gf-hdg-dev` (trunk) | none |
| `gf-hdg-subdomains-dev` | none |
| `gf-hdg-p-adaptivity` | none |
| `gf-hdg-linearise-first` | **yes** |
| `gf-interp-hdg-dev` | **yes** |

The trunk's `ReconstructTotalFlux` is the scalar one and aliases `b_z` itself,
which is registered and safe. The seven views arrived with the `neq > 1`
generalisation. **meq is affected**, because `meq-integration` merges lf — so
the practical answer to you is unchanged — but a reader of your file would
otherwise go looking for it on branches that never had it.

## And your "why nobody has met this" is right by the wrong route

You attribute `cuda`'s immunity to `mt_host` being true so the destructor is
inert. Since `h_mt` decides nothing, that cannot be it. The eviction runs under
`cuda` too — it just erases a key that is not there, because with
`GetHostMemoryType() == HOST` the **owner** was never registered either. Same
conclusion, and worth having right, because it says the immunity is a property
of the *host* memory type and not of the delete path.

## What is left open, and it is upstream's rather than ours

Two questions this fix does not answer, and we have not taken either:

* should `Register_` refuse, or at least warn, on an address already in
  `maps->memories`? `Insert`'s `emplace` makes the collision silent outside
  `MFEM_DEBUG`, and that silence is what turns a double registration into a
  later abort two layers away;
* should a **lazily** registered `Memory` — one registered by `MakeAlias` or by
  a device-class `Read`/`Write` rather than by its own constructor — claim
  `OWNS_INTERNAL`? It is the `OWNS_INTERNAL` that licenses the `Erase`.

Either would close the trap for the next caller, which is the thing your third
option correctly says a local fix does not do. Both are changes to `Memory`'s
ownership convention in MFEM core, which is a separate branch off master with
its own justification — the shape `pardiso-device-host-sync` had — rather than
something to smuggle in through `fem/darcy`. If you would rather it were closed
upstream than worked around per caller, say so and it goes on the list with your
name on the reason.

## Thanks for the instrument

`Device("debug")` plus a breakpoint on the protector is now twice the thing that
turned a day into an afternoon here, and it would have been the third if I had
reached for it instead of reading `Delete_` and stopping. Your report had the
abort, the register values and the two-hit count in it; everything I had to add
was two arms of a probe.
