// external crates
extern crate postgres;

use std::collections::HashMap;
use std::hash::Hash;

use proteomic::models::persistable::Persistable;

pub trait Collectable {
    fn get_collection_identifier(&self) -> &String;
}

pub struct Collection<V> where V: PartialEq + Eq + Hash + Persistable<V, i32, String> + Collectable {
    objects: HashMap<String, V>
}

impl<V> Collection<V> where V: PartialEq + Eq + Hash + Persistable<V, i32, String> + Collectable{

    pub fn new() -> Collection<V> where V: PartialEq + Eq + Hash + Persistable<V, i32, String> + Collectable {
        return Collection {
            objects: HashMap::new()
        }
    }

    pub fn add(&mut self, item: V) -> bool {
        if !self.contains(&item) {
            self.objects.insert(item.get_collection_identifier().clone(), item);
            return true;
        }
        return false;
    }

    pub fn remove(&mut self, item: &V) {
        self.objects.remove(item.get_collection_identifier());
    }

    pub fn clear(&mut self) {
        self.objects.clear();
    }

    pub fn contains(&self, item: &V) -> bool {
        return self.objects.contains_key(item.get_collection_identifier());
    }

    pub fn len(&self) -> usize {
        return self.objects.len();
    }

    pub fn get(&self, identifier: &String) -> Option<&V>{
        return self.objects.get(identifier);
    }

    // pub fn iter_mut(&self) -> IterMut<String, V> {
    //     self.objects.iter_mut();
    // }

    // pub fn save(&self, conn: &postgres::Connection) {
    //     let transaction = conn.transaction().unwrap();
    //     let prepared_insert_statement = transaction.prepare_cached(T::get_insert_query()).unwrap();
    //     //let prepared_update_statement = transaction.prepare_cached(T::get_update_statement()).unwrap();
    //     for (_key, t_object) in self.objects.iter() {
    //         // if t_object.get_id() > 0 {
    //         //     t_object.execute_update_statement(&prepared_update_statement);
    //         // } else {
    //         //     t_object.execute_insert_statement(&prepared_insert_statement);
    //         // }
    //         t_object.execute_insert_statement(&prepared_insert_statement);
    //     }
    //     transaction.commit();
    //     transaction.finish();
    // }
}