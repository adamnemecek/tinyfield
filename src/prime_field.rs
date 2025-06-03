use core::{cmp, convert, fmt, marker, ops};

/// A representation of a field with prime order.
///
/// This trait is implemented by any structure purporting to be a prime
/// field.  Said structures, such as `GF2` or `GF3` are exposed in the `field`
/// module and applications should not implement their own.
pub trait PrimeField: marker::Sized + core::fmt::Debug + marker::Copy {
    /// The characteristic of this field, a prime.
    ///
    /// The characteristic of an _element_ of a ring is the fewest times it
    /// must be added to itself in order to reach zero.  In a field, this
    /// number is the same for each element of the structure.  Further it is a
    /// prime number.  Thus one may talk of the characteristic of a field.
    ///
    /// This is the only information required to uniquely determine this field.
    /// It is also the size of this field.
    const CHARACTERISTIC: u8;

    /// A map of multiplicative inverses for this field.
    ///
    // TODO(robert) Evaluate replacing this with Euclid's algorithm for speed.
    const DIVISION_TABLE: [u8; 256];

    /// The zero element of the field.
    const zero: PrimeFieldElt<Self> = PrimeFieldElt::new(0);

    /// The multiplicative unit of the field.
    const one: PrimeFieldElt<Self> = PrimeFieldElt::new(1);

    /// An iterator over the p elements of this field.
    ///
    // TODO(robert) Construct a specific `PrimeFieldIterator` or investigate
    //              implementing iter on GFp
    fn elts() -> core::iter::Scan<
        ops::Range<u8>,
        PrimeFieldElt<Self>,
        fn(&mut PrimeFieldElt<Self>, u8) -> Option<PrimeFieldElt<Self>>,
    > {
        (0..Self::CHARACTERISTIC).scan(Self::zero, |acc, _| Some(*acc + Self::one))
    }
}

/// An element of a prime field.
///
/// Prime field elements should never be instanciated by external code but
/// instead either created by iterating over a field or using the built-in
/// constants `PrimeField::zero` and `PrimeField::one`.
///
/// Represented internally as a single byte.
#[derive(Clone, Copy)]
pub struct PrimeFieldElt<F: PrimeField> {
    val: u8,
    phantom: marker::PhantomData<F>,
}

impl<F: PrimeField> PrimeFieldElt<F> {
    #[inline]
    pub const fn new(val: u8) -> Self {
        Self {
            val,
            phantom: marker::PhantomData,
        }
    }
}

impl<F: PrimeField> fmt::Debug for PrimeFieldElt<F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_fmt(format_args!("F{}({})", F::CHARACTERISTIC, self.val))
    }
}

impl<F: PrimeField> ops::Add for PrimeFieldElt<F> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        Self::new((((self.val as u16) + (rhs.val as u16)) % (F::CHARACTERISTIC as u16)) as u8)
    }
}

impl<F: PrimeField> ops::Neg for PrimeFieldElt<F> {
    type Output = Self;

    fn neg(self) -> Self {
        Self::new((F::CHARACTERISTIC - self.val) % F::CHARACTERISTIC)
    }
}

impl<F: PrimeField> ops::Sub for PrimeFieldElt<F> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        self + (-rhs)
    }
}

impl<F: PrimeField> ops::Mul for PrimeFieldElt<F> {
    type Output = Self;

    fn mul(self, rhs: PrimeFieldElt<F>) -> PrimeFieldElt<F> {
        &self * &rhs
    }
}

impl<F: PrimeField> ops::Mul for &PrimeFieldElt<F> {
    type Output = PrimeFieldElt<F>;

    fn mul(self, rhs: Self) -> Self::Output {
        Self::Output::new(
            (((self.val as u16) * (rhs.val as u16)) % (F::CHARACTERISTIC as u16)) as u8,
        )
    }
}

impl<F: PrimeField> ops::Div for PrimeFieldElt<F> {
    type Output = Self;

    fn div(self, rhs: Self) -> Self {
        assert_ne!(rhs, F::zero, "Division by zero");
        self * Self::new(F::DIVISION_TABLE[rhs.val as usize])
    }
}

impl<F: PrimeField> PrimeFieldElt<F> {
    pub fn pow(self, rhs: u8) -> Self {
        let rhs = rhs % (F::CHARACTERISTIC - 1);
        if rhs == 0 {
            F::one
        } else {
            self * self.pow(rhs - 1)
        }
    }
}

impl<F: PrimeField> cmp::PartialEq for PrimeFieldElt<F> {
    fn eq(&self, other: &Self) -> bool {
        self.val == other.val
    }
}

impl<F: PrimeField> cmp::Eq for PrimeFieldElt<F> {}

impl<F: PrimeField> convert::From<u8> for PrimeFieldElt<F> {
    fn from(x: u8) -> Self {
        Self {
            val: x % F::CHARACTERISTIC,
            phantom: marker::PhantomData,
        }
    }
}

impl<F: PrimeField> convert::From<PrimeFieldElt<F>> for u8 {
    fn from(x: PrimeFieldElt<F>) -> Self {
        x.val
    }
}
